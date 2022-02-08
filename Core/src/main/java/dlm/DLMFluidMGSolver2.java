package dlm;

import basic.ScalarFunction;
import basic.VectorFunctionOnCells;
import io.vavr.Tuple2;
import it.unimi.dsi.fastutil.ints.Int2DoubleMap;
import linalg.*;
import mixed.*;
import multigrid.MGPreconditionerSpace;
import multigrid.Smoother;
import org.jetbrains.annotations.NotNull;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class DLMFluidMGSolver2
	extends DLMSolver
{
	public double smootherOmega = 1;
	public int smootherRuns = 3;
	
	RichardsonExplicitSchur schur;
	MultiGridFluid fluid;
	
	public DLMFluidMGSolver2(final MultiGridFluid f)
	{
		fluid = f;
	}
	
	@NotNull
	private MGPreconditionerSpace<TaylorHoodSpace, TPCell, TPFace, QkQkFunction, MixedValue, MixedGradient, MixedHessian> create_space(
		final double dt,
		final double t,
		final FluidIterate iterate)
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
				velocity = new MixedTPFESpaceFunction<>(space.getShapeFunctionMap(),
				                                        restrictedIterate).getVelocityFunction();
				final FluidSystem fs = fluid.getFluidSystemForSpace(space, velocity, 0,
				                                                    restrictedIterate);
				final var blockRhs = Fluid.getBlockRhs(fs, dt);
				final SparseMatrix s = new SparseMatrix(blockRhs._1);
				final Int2DoubleMap nodeValues = fluid.getDirichletNodeValuesForSpace(space, t);
				nodeValues.forEach((node, val) ->
				                   {
					                   s.deleteColumn(node);
					                   s.deleteRow(node);
					                   s.set(1, node, node);
				                   });
				return new Tuple2<>(s, new DenseVector(n));
			}
			
			@Override
			public List<Smoother> createSmoothers()
			{
				verbose = true;
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
				                                fluid.getDirichletBoundary(),
				                                (f, fun) -> fun.hasVelocityFunction(),
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
			schur = new RichardsonExplicitSchur(systemMatrix);
		else
		{
			System.out.println("reset");
			schur.resetOffDiagonals(systemMatrix);
		}
		final MGPreconditionerSpace<TaylorHoodSpace, TPCell, TPFace, QkQkFunction, MixedValue, MixedGradient, MixedHessian>
			mg = create_space(dt, t, fluidState);
		final Smoother Sgs = new UzawaStokesSmoother3(10,
		                                              1,
		                                              mg.getFinestSpace()
		                                                .getVelocitySize());
		mg.verbose = false;
		schur.preconditioner = new VectorMultiplyable()
		{
			@Override
			public int getVectorSize()
			{
				return mg.getVectorSize();
			}
			
			@Override
			public int getTVectorSize()
			{
				return mg.getTVectorSize();
			}
			
			@Override
			public Vector mvMul(final Vector vector)
			{
				DenseVector rhs = new DenseVector(vector);
				Vector mgSol = vector.mul(0);
				Vector force = new DenseVector(vector);
				for (int j = 0; j < 1; j++)
				{
					final Vector correction = mg.mvMul(rhs);
					mgSol = mgSol.add(correction);
					System.out.println("MGSOL " + mg.finest_system
						.mvMul(mgSol)
						.sub(force)
						.euclidianNorm());
					System.out.println("MGSOL " + mg.finest_system
						.mvMul(correction)
						.sub(rhs)
						.euclidianNorm());
					System.out.println("MGSOL " + schur.getSchurComplement()
					                                   .mvMul(mgSol)
					                                   .sub(vector)
					                                   .euclidianNorm());
					mgSol = Sgs.smooth(schur.getSchurComplement(),
					                   vector,
					                   mgSol,
					                   true,
					                   "SGS   ");
					System.out.println("SGSOL " + schur.getSchurComplement()
					                                   .mvMul(mgSol)
					                                   .sub(vector)
					                                   .euclidianNorm());
					final Vector finalFluid = mgSol;
					final List<Vector> forces
						= IntStream.range(0, particleStates.size())
						           .mapToObj(i -> schur.getTopBlock(i)
						                               .mvMul(schur.getDiagonalInverse(i)
						                                           .mvMul(schur.getLeftBlock(i)
						                                                       .mvMul(finalFluid))))
						           .collect(Collectors.toList());
					force = new DenseVector(vector);
					System.out.println(force.euclidianNorm());
					for (final Vector v : forces)
						force = force.add(v);
					System.out.println("FFF " + mg.finest_system.mvMul(mgSol)
					                                            .sub(force)
					                                            .euclidianNorm());
					System.out.println(force.euclidianNorm());
					rhs = new DenseVector(force.sub(mg.finest_system.mvMul(mgSol)));
					System.out.println("RESIUDIAL " + rhs.euclidianNorm());
				}
				return mgSol;
//				Sgs.smooth(schur.getSchurComplement(),
//				                  vector,
//				                  mgSol,
//				                  true,
//				                  "SGS   ");
			}
			
			@Override
			public Vector tvMul(final Vector vector)
			{
				throw new UnsupportedOperationException("not implemented yet");
			}
		};
		return schur.mvMul(rhs);
	}
}
