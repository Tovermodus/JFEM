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

public class DLMSplittingSolver
	extends DLMSolver
{
	public double splittingFactor = 1;
	RichardsonExplicitSchur schur;
	MultiGridFluid fluid;
	
	public DLMSplittingSolver(final MultiGridFluid f)
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
		mg.verbose = true;
		mg.cycles = 8;
		final Smoother Sgs = new BSSmoother2(1,
		                                     1.0,
		                                     mg.getFinestSpace()
		                                       .getVelocitySize());
		final VectorMultiplyable particleForces = new VectorMultiplyable()
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
				final List<Vector> forces
					= IntStream.range(0, particleStates.size())
					           .mapToObj(i -> schur.getTopBlock(i)
					                               .mvMul(schur.getDiagonalInverse(i)
					                                           .mvMul(schur.getLeftBlock(i)
					                                                       .mvMul(vector))))
					           .collect(Collectors.toList());
				Vector force = new DenseVector(vector.getLength());
				for (final Vector v : forces)
					force = force.add(v);
				force = force.add(mg.finest_system.mvMul(vector)
				                                  .mul(splittingFactor - 1));
				return force;
			}
			
			@Override
			public Vector tvMul(final Vector vector)
			{
				final List<Vector> forces
					= IntStream.range(0, particleStates.size())
					           .mapToObj(i -> schur.getTopBlock(i)
					                               .mvMul(schur.getDiagonalInverse(i)
					                                           .tvMul(schur.getLeftBlock(i)
					                                                       .mvMul(vector))))
					           .collect(Collectors.toList());
				Vector force = new DenseVector(vector.getLength());
				for (final Vector v : forces)
					force = force.add(v);
				force = force.add(mg.finest_system.tvMul(vector)
				                                  .mul(splittingFactor - 1));
				return force;
			}
		};
		System.out.println();
		System.out.println();
		System.out.println();
		System.out.println();
		System.out.println();
		final SparseMatrix F = ((SparseMatrix) mg.finest_system).mul(splittingFactor);
		final double cbceig = particleForces.powerIterationNonSymmetric();
		final double finveig = F.inverse()
		                        .powerIterationNonSymmetric();
		System.out.println("#####     CBC " + cbceig);
		System.out.println("#####     Finv " + finveig);
		System.out.println("#####     F " + F.powerIterationNonSymmetric());
		System.out.println("#####     A " + F.addVm(particleForces.mulVm(-1))
		                                     .powerIterationNonSymmetric());
		System.out.println("#####     A " + schur.getSchurComplement()
		                                         .powerIterationNonSymmetric());
		System.out.println("#####     FinvCBC " + VectorMultiplyable.concatenate(F.inverse(),
		                                                                         particleForces)
		                                                            .powerIterationNonSymmetric());

//		for (int i = 0; i < schur.getBlockSize() - 1; i++)
//		{
//			final int size = schur.getDiagonalInverse(i)
//			                      .getVectorSize() / 2;
//			Matrix blockinv = new DenseMatrix(schur.getDiagonalBlock(i)
//			                                       .slice(new IntCoordinates(0, 0),
//			                                              new IntCoordinates(size, size))
//			                                       .inverse());
//			System.out.println("#####    Minv " + blockinv
//				.powerIterationNonSymmetric());
//			blockinv = new DenseMatrix(schur.getDiagonalBlock(i)
//			                                .slice(new IntCoordinates(0, 0),
//			                                       new IntCoordinates(size, size)));
//			System.out.println("#####    M " + blockinv
//				.powerIterationNonSymmetric());
//			blockinv = new DenseMatrix(schur.getDiagonalBlock(i)
//			                                .slice(new IntCoordinates(0, size),
//			                                       new IntCoordinates(size, 2 * size))
//			                                .inverse());
//			System.out.println("#####    Cinv " + blockinv.powerIterationNonSymmetric());
//			blockinv = new DenseMatrix(schur.getDiagonalBlock(i)
//			                                .inverse());
//			System.out.println("#####    Blockinv" + blockinv
//				.powerIterationNonSymmetric());
//			blockinv = new DenseMatrix(schur.getDiagonalBlock(i)
//			                                .inverse()
//			                                .slice(new IntCoordinates(size, size),
//			                                       new IntCoordinates(size * 2, size * 2)));
//			System.out.println("#####    Sinv " + blockinv
//				.powerIterationNonSymmetric());
//			System.out.println();
//		}
//		final SparseMatrix ADiag = new SparseMatrix(new BlockDenseMatrix((SparseMatrix) mg.finest_system,
//		                                                                 mg.getVectorSize()).getInvertedDiagonalMatrix());
//		System.out.println("##########################   A " + ((SparseMatrix) mg.finest_system).mul(
//			                                                                                        splittingFactor)
//		                                                                                        .powerIterationNonSymmetric());
//		System.out.println("##########################   Ainv " + ((SparseMatrix) mg.finest_system).mul(
//			                                                                                           splittingFactor)
//		                                                                                           .inverse()
//		                                                                                           .powerIterationNonSymmetric());
//		System.out.println("##########################   K " + particleForces.powerIterationNonSymmetric());
//		System.out.println("##########################   AinvK " + VectorMultiplyable.concatenate(((SparseMatrix) mg.finest_system).mul(
//			                                                                                                                           splittingFactor)
//		                                                                                                                           .inverse(),
//		                                                                                          particleForces)
//		                                                                             .powerIterationNonSymmetric());
		System.out.println();
		System.out.println();
		System.out.println();
		System.out.println();
		System.out.println();
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
				Vector rhs = new DenseVector(vector);
				Vector iterate = vector.mul(0);
				for (int j = 0; j < 5; j++)
				{
					System.out.println();
					System.out.println();
					final Vector newRhs = vector.add(particleForces.mvMul(iterate));
					System.out.println("rhs difference " + newRhs.sub(rhs)
					                                             .euclidianNorm() + "    " +
						                   "  b+CBCx2 - (b+CBCx0)");
					System.out.println("max iterate difference " + finveig * newRhs.sub(rhs)
					                                                               .euclidianNorm()
						                   + "      x1-x0");
					rhs = newRhs;
					System.out.println("mg presolution residual" + F.mvMul(iterate)
					                                                .sub(rhs)
					                                                .euclidianNorm()
						                   + "       Fx0 = b+CBCx0");
					mg.verbose = false;
					final Vector newiterate = mg.mvMul(rhs)
					                            .mul(1. / splittingFactor);
					System.out.println("mg solution residual" + F.mvMul(newiterate)
					                                             .sub(rhs)
					                                             .euclidianNorm() + "      Fx1 = b+CBCx0");
					System.out.println("iterate difference " + newiterate.sub(iterate)
					                                                     .euclidianNorm() + "      x1-x0");
					System.out.println("max residual difference " + cbceig * newiterate.sub(iterate)
					                                                                   .euclidianNorm());
					iterate = newiterate;
					System.out.println("schur residual " + schur.getSchurComplement()
					                                            .mvMul(iterate)
					                                            .sub(vector)
					                                            .euclidianNorm()
						                   + "      (F-CBC)x1 = b");
					iterate = Sgs.smooth(schur.getSchurComplement(), vector, iterate, true,
					                     "              SPLITTING ");
					System.out.println("schur residual " + schur.getSchurComplement()
					                                            .mvMul(iterate)
					                                            .sub(vector)
					                                            .euclidianNorm()
						                   + "      (F-CBC)x2 = b");
				}
				return iterate;
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
