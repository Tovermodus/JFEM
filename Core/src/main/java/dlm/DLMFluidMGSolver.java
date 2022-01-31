package dlm;

import basic.ScalarFunction;
import basic.VectorFunction;
import basic.VectorFunctionOnCells;
import io.vavr.Tuple2;
import linalg.*;
import mixed.*;
import multigrid.MGPreconditionerSpace;
import multigrid.Smoother;
import org.jetbrains.annotations.NotNull;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

import java.util.ArrayList;
import java.util.List;

public class DLMFluidMGSolver
	extends DLMSolver
{
	public double smootherOmega = 1;
	public int smootherRuns = 7;
	
	PreconditionedIterativeImplicitSchur schur;
	MultiGridFluid fluid;
	
	public DLMFluidMGSolver(final MultiGridFluid f)
	{
		fluid = f;
	}
	
	@NotNull
	private MGPreconditionerSpace<TaylorHoodSpace, TPCell, TPFace, QkQkFunction, MixedValue, MixedGradient, MixedHessian> create_space(
		final MultiGridFluid fluid, final double dt, final double t, final FluidIterate iterate)
	{
		final MGPreconditionerSpace<TaylorHoodSpace, TPCell, TPFace, QkQkFunction, MixedValue, MixedGradient, MixedHessian>
			mg = new MGPreconditionerSpace<>(fluid.refinements, fluid.polynomialDegree)
		{
			@Override
			public List<TaylorHoodSpace> createSpaces(final int refinements)
			{
				final List<TaylorHoodSpace> ret = new ArrayList<>();
				for (int i = 0; i < refinements + 1; i++)
				{
					final TaylorHoodSpace subSpace = new TaylorHoodSpace(fluid.startCoordinates,
					                                                     fluid.endCoordinates,
					                                                     fluid.coarsestCells.mul((int) Math.pow(
						                                                     2,
						                                                     i)));
					subSpace.assembleCells();
					subSpace.assembleFunctions(fluid.polynomialDegree);
					ret.add(subSpace);
				}
				return ret;
			}
			
			@Override
			public Tuple2<VectorMultiplyable, DenseVector> createSystem(final TaylorHoodSpace space)
			{
				final int n = space.getShapeFunctionMap()
				                   .size();
				final VectorFunctionOnCells<TPCell, TPFace> velocity;
				final Vector restrictedIterate;
				if (iterate != null)
				{
					restrictedIterate = restrictToSize(n, iterate.current);
				} else
				{
					restrictedIterate = new DenseVector(n);
				}
				velocity = new MixedTPFESpaceFunction<QkQkFunction>(space.getShapeFunctionMap(),
				                                                    restrictedIterate).getVelocityFunction();
				final FluidSystem fs = fluid.getFluidSystemForSpace(space, velocity, 0,
				                                                    restrictedIterate);
				final int firstPRessure = (int) space.getShapeFunctionMap()
				                                     .values()
				                                     .stream()
				                                     .filter(ComposeMixedShapeFunction::hasVelocityFunction)
				                                     .count();
				final var blockRhs = fluid.getBlockRhsForSpace(space, firstPRessure, fs, dt, t);
				final SparseMatrix s = new SparseMatrix(blockRhs._1);
				space.writeBoundaryValuesTo(new ComposedMixedFunction(VectorFunction.fromLambda(
					                            fluid.velocityBoundaryValues(t),
					                            2,
					                            2)), s,
				                            new DenseVector(n));
				space.overWriteValue(firstPRessure, 0, s, new DenseVector(n));
				System.out.println("diffffference" + s.sub(blockRhs._1)
				                                      .absMaxElement());
				return new Tuple2<>(blockRhs._1, blockRhs._2);
			}
			
			@Override
			public List<Smoother> createSmoothers()
			{
				verbose = false;
				final List<Smoother> ret = new ArrayList<>();
				for (int i = 1; i < fluid.refinements + 1; i++)
				{
					ret.add(new BSSmoother2(smootherRuns,
					                        smootherOmega,
					                        spaces.get(i)
					                              .getVelocitySize()));
				}
				return ret;
			}
			
			@Override
			public void applyZeroBoundaryConditions(final TaylorHoodSpace space, final MutableVector vector)
			{
				space.projectOntoBoundaryValues(new ComposedMixedFunction(ScalarFunction.constantFunction(
					                                                                        0)
				                                                                        .makeIsotropicVectorFunction()),
				                                vector);
				vector.set(0, space.getVelocitySize());
			}
		};
		final DenseVector v1 = new DenseVector(mg.systems.get(0)
		                                                 .getVectorSize());
		final DenseVector v2 = new DenseVector(mg.systems.get(0)
		                                                 .getVectorSize());
		final VectorMultiplyable P = mg.prolongationMatrices.get(0);
		int j = 6;
		for (final IntCoordinates c : v1.getShape()
		                                .range())
		{
			v1.set(j++, c);
			v2.set(1, c);
		}
		mg.applyZeroBoundaryConditions(mg.spaces.get(0), v1);
		mg.applyZeroBoundaryConditions(mg.spaces.get(0), v2);
		System.out.println("subspacesolver   " + Math.abs(mg.systems.get(1)
		                                                            .mvMul(P.mvMul(v1))
		                                                            .inner(P.mvMul(v2)) - mg.systems.get(0)
		                                                                                            .mvMul(v1)
		                                                                                            .inner(v2)));
		return mg;
	}
	
	@Override
	protected Vector solve(final BlockSparseMatrix systemMatrix,
	                       final DenseVector rhs,
	                       final FluidIterate fluidState,
	                       final List<ParticleIterate> particleStates,
	                       final FluidSystem fluidSystem,
	                       final List<ParticleSystem> particleSystems, final double dt, final double t)
	{
		
		if (schur == null)
			schur = new PreconditionedIterativeImplicitSchur(systemMatrix,
			                                                 create_space(fluid, dt, t, fluidState));
		else
			schur.resetOffDiagonals(systemMatrix);
		return schur.mvMul(rhs);
	}
}
