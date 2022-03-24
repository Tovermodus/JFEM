package dlm;

import basic.ScalarFunction;
import basic.VectorFunctionOnCells;
import io.vavr.Tuple2;
import it.unimi.dsi.fastutil.ints.Int2DoubleMap;
import linalg.*;
import mixed.*;
import multigrid.BlockGaussseidelSmoother;
import multigrid.MGPreconditionerSpace;
import multigrid.Smoother;
import org.jetbrains.annotations.NotNull;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ConcurrentHashMap;
import java.util.stream.IntStream;

public class DLMDiagonal
	extends DLMSolver
{
	MultiGridFluid fluid;
	
	public DLMDiagonal(final MultiGridFluid f)
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
				final var blockRhs = Fluid.assembleBlockRhs(fs, dt);
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
					ret.add(new BSSmoother2(6,
					                        1,
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
	
	ConcurrentHashMap<Integer, DenseMatrix> inverses;
	
	DenseMatrix getDiagonalInverse(final BlockSparseMatrix mat, final int i)
	{
		synchronized (this)
		{
			if (inverses == null)
			{
				inverses = new ConcurrentHashMap<>();
				IntStream.range(0, mat.getBlockEnds().length - 1)
				         .parallel()
				         .forEach(j -> inverses.put(j,
				                                    mat.getBlockMatrix(j + 1, j + 1)
				                                       .inverse()));
			}
		}
		return inverses.get(i);
	}
	
	@Override
	protected Vector solve(final BlockSparseMatrix systemMatrix,
	                       final DenseVector rhs,
	                       final FluidIterate fluidState,
	                       final List<ParticleIterate> particleStates,
	                       final FluidSystem fluidSystem,
	                       final List<ParticleSystem> particleSystems, final double dt, final double t)
	{
		
		final int[] blockStarts = systemMatrix.getBlockStarts();
		final int[] blockEnds = systemMatrix.getBlockEnds();
		final MGPreconditionerSpace<TaylorHoodSpace, TPCell, TPFace, QkQkFunction, MixedValue, MixedGradient, MixedHessian>
			mg = create_space(dt, t, fluidState);
		mg.verbose = true;
		final VectorMultiplyable precond = new VectorMultiplyable()
		{
			@Override
			public int getVectorSize()
			{
				return systemMatrix.getVectorSize();
			}
			
			@Override
			public int getTVectorSize()
			{
				return systemMatrix.getTVectorSize();
			}
			
			@Override
			public Vector mvMul(final Vector vector)
			{
//				final List<DenseVector> subVectors = BlockSparseMatrix.partitionVector(vector,
//				                                                                       blockStarts,
//				                                                                       blockEnds);
//				final List<DenseVector> preconditionedVectors =
//					IntStream.range(1,
//					                particleStates.size())
//					         .parallel()
//					         .mapToObj(i ->
//					                   {
//						                   return getDiagonalInverse(systemMatrix,
//						                                             i - 1).mvMul(
//							                   subVectors.get(i));
//					                   })
//					         .collect(Collectors.toList());
//				final DenseVector fluidRhs = new DenseVector(subVectors.get(0));
//				IntStream.range(1,
//				                particleStates.size())
//				         .forEach(i -> fluidRhs.subInPlace(systemMatrix.getBlockMatrix(0, i)
//				                                                       .mvMul(preconditionedVectors.get(
//					                                                       i - 1))));
//				final DenseVector ret = new DenseVector(vector.getLength());
//				ret.addSmallVectorAt(systemMatrix.getBlockMatrix(0, 0)
//				                                 .solve(fluidRhs), 0);
//				for (int i = 1; i < particleStates.size(); i++)
//					ret.addSmallVectorAt(preconditionedVectors.get(i - 1), blockStarts[i]);
//
//				final List<DenseVector> residualSubVectors =
//					BlockSparseMatrix.partitionVector(vector.sub(systemMatrix.mvMul(ret)),
//					                                  blockStarts,
//					                                  blockEnds);
//				System.out.println();
//
//				residualSubVectors.forEach(v -> System.out.println(v.euclidianNorm()));
//				System.out.println(vector.sub(systemMatrix.mvMul(ret))
//				                         .euclidianNorm() + " " + vector.euclidianNorm());
//				System.out.println();
//				return ret;
				final Vector ret = new BlockGaussseidelSmoother(systemMatrix, 3).smooth(systemMatrix,
				                                                                        vector,
				                                                                        vector,
				                                                                        true,
				                                                                        "GS ");
				System.out.println(systemMatrix.mvMul(ret)
				                               .sub(vector)
				                               .euclidianNorm());
				return ret;
			}
			
			@Override
			public Vector tvMul(final Vector vector)
			{
				throw new UnsupportedOperationException("not implemented yet");
			}
		};
		final IterativeSolver it = new IterativeSolver(true);
		it.showProgress = true;
		return it.solvePGMRES(systemMatrix, precond, rhs, 1e-6);
	}
}
