package dlm;

import basic.ScalarFunction;
import basic.VectorFunction;
import basic.VectorFunctionOnCells;
import com.google.common.base.Stopwatch;
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

public abstract class MultiGridFluid
	implements Fluid
{
	public MGPreconditionerSpace<TaylorHoodSpace, TPCell, TPFace, QkQkFunction, MixedValue, MixedGradient,
		MixedHessian>
		space;
	private final CoordinateVector startCoordinates;
	private final CoordinateVector endCoordinates;
	private final IntCoordinates coarsestCells;
	private final int refinements;
	private final int polynomialDegree;
	private final double dt;
	
	public MultiGridFluid(final CoordinateVector startCoordinates, final CoordinateVector endCoordinates,
	                      final IntCoordinates coarsestCells, final int refinements, final int polynomialDegree,
	                      final double dt)
	{
		this.startCoordinates = startCoordinates;
		this.endCoordinates = endCoordinates;
		this.coarsestCells = coarsestCells;
		this.refinements = refinements;
		this.polynomialDegree = polynomialDegree;
		this.dt = dt;
		space = create_space(startCoordinates,
		                     endCoordinates,
		                     coarsestCells,
		                     refinements,
		                     polynomialDegree,
		                     dt, null);
	}
	
	@NotNull
	private MGPreconditionerSpace<TaylorHoodSpace, TPCell, TPFace, QkQkFunction, MixedValue, MixedGradient, MixedHessian> create_space(
		final CoordinateVector startCoordinates,
		final CoordinateVector endCoordinates,
		final IntCoordinates coarsestCells,
		final int refinements,
		final int polynomialDegree,
		final double dt, final FluidIterate iterate)
	{
		final MGPreconditionerSpace<TaylorHoodSpace, TPCell, TPFace, QkQkFunction, MixedValue, MixedGradient, MixedHessian>
			mg = new MGPreconditionerSpace<>(refinements, polynomialDegree)
		{
			@Override
			public List<TaylorHoodSpace> createSpaces(final int refinements)
			{
				final List<TaylorHoodSpace> ret = new ArrayList<>();
				for (int i = 0; i < refinements + 1; i++)
				{
					final TaylorHoodSpace subSpace = new TaylorHoodSpace(startCoordinates,
					                                                     endCoordinates,
					                                                     coarsestCells.mul((int) Math.pow(
						                                                     2,
						                                                     i)));
					subSpace.assembleCells();
					subSpace.assembleFunctions(polynomialDegree);
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
				final FluidSystem fs = getFluidSystemForSpace(space, velocity, 0, restrictedIterate);
				final int firstPRessure = (int) space.getShapeFunctionMap()
				                                     .values()
				                                     .stream()
				                                     .filter(ComposeMixedShapeFunction::hasVelocityFunction)
				                                     .count();
				final var blockRhs = getBlockRhsForSpace(space, firstPRessure, fs, dt);
				final SparseMatrix s = new SparseMatrix(blockRhs._1);
				space.writeBoundaryValuesTo(new ComposedMixedFunction(VectorFunction.fromLambda(
					                            velocityBoundaryValues(),
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
				for (int i = 1; i < refinements + 1; i++)
				{
//					ret.add(new UzawaStokesSmoother2(5,
//					                                 0.5,
//					                                 (int) spaces.get(i)
//					                                             .getShapeFunctionMap()
//					                                             .values()
//					                                             .stream()
//					                                             .filter(ComposeMixedShapeFunction::hasVelocityFunction)
//					                                             .count()));
					//ret.add(new RichardsonSmoother(0.05, 4));
					ret.add(new BSSmoother2(3,
					                        0.5,
					                        spaces.get(i)
					                              .getVelocitySize()));
				}
				return ret;
			}
			
			@Override
			public void applyCorrectBoundaryConditions(final TaylorHoodSpace space,
			                                           final MutableVector vector)
			{
				//throw new UnsupportedOperationException("not implemented yet");
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
	public TaylorHoodSpace getSpace()
	{
		return space.spaces.get(refinements);
	}
	
	@Override
	public FluidSystem buildSystem(final double t, final FluidIterate iterate)
	{
		final FluidSystem ret = getFluidSystemForSpace(getSpace(), getVelocity(iterate), t, iterate.current);
		space = create_space(startCoordinates,
		                     endCoordinates,
		                     coarsestCells,
		                     refinements,
		                     polynomialDegree,
		                     dt, iterate);
		return ret;
	}
	
	@NotNull
	private FluidSystem getFluidSystemForSpace(final TaylorHoodSpace space, final VectorFunctionOnCells<TPCell,
		TPFace> velocity, final double t,
	                                           final Vector currentIterate)
	{
		final int n = space.getShapeFunctionMap()
		                   .size();
		final SparseMatrix massMatrix = new SparseMatrix(n, n);
		space.writeCellIntegralsToMatrix(getMassIntegrals(), massMatrix);
		final SparseMatrix flowMatrix = new SparseMatrix(n, n);
		space.writeCellIntegralsToMatrix(getIntegrals(), flowMatrix);
		final SparseMatrix semiImplicitMatrix = new SparseMatrix(n, n);
		System.out.println("semiImp");
		final Stopwatch s = Stopwatch.createStarted();
		space.writeCellIntegralsToMatrix(getSemiImplicitIntegrals(velocity),
		                                 semiImplicitMatrix);
		System.out.println("semiImpDone" + s.elapsed());
		final DenseVector forceRhs = new DenseVector(n);
		space.writeCellIntegralsToRhs(getForceIntegrals(t), forceRhs);
		final DenseVector accelRhs = massMatrix.mvMul(currentIterate);
		return new FluidSystem(massMatrix, flowMatrix, semiImplicitMatrix, forceRhs, accelRhs);
	}
}
