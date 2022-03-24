package dlm;

import basic.VectorFunctionOnCells;
import io.vavr.Tuple2;
import it.unimi.dsi.fastutil.ints.Int2DoubleMap;
import linalg.*;
import mixed.*;
import multigrid.AMGPreconditionerSpace;
import multigrid.ForwardBackwardGaussSeidelSmoother;
import multigrid.Smoother;
import org.jetbrains.annotations.NotNull;
import tensorproduct.ContinuousTPFEVectorSpace;
import tensorproduct.ContinuousTPVectorFunction;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class DLMLagrangeAMGSolver
	extends DLMSolver
{
	RichardsonSchur schur;
	MultiGridFluid fluid;
	final List<Particle> particles;
	
	public DLMLagrangeAMGSolver(final MultiGridFluid f, final List<Particle> particles)
	{
		fluid = f;
		this.particles = particles;
	}
	
	@NotNull
	private DLMAMG create_space(
		final List<ParticleIterate> particleStates, final double dt, final double t,
		final FluidIterate iterate)
	{
		return new DLMAMG(fluid.refinements, fluid.polynomialDegree)
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
			public Tuple2<SparseMatrix, DenseVector> createSystem(final TaylorHoodSpace space)
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
				final List<Matrix> schurContributions =
					IntStream.range(0, particles.size())
					         .parallel()
					         .mapToObj(i -> new Tuple2<>(i,
					                                     particles.get(i)
					                                              .buildLagrangeBackgroundMatrix(
						                                              space,
						                                              particleStates.get(i))))
					         .map(e ->
					              {
						              final SparseMatrix block =
							              new SparseMatrix(s.getRows(),
							                               2 * e._2.getCols());
						              block.addSmallMatrixInPlaceAt(e._2,
						                                            0,
						                                            e._2.getCols());
						              return new Tuple2<>(e._1, block);
					              })
					         .map(e ->
						              e._2.mmMul(new SparseMatrix(schur.getDiagonalInverse(
							               e._1)))
						                  .mtMul(e._2)
						                  .mul(-1))
					         .collect(
						         Collectors.toList());
				for (final Matrix m : schurContributions)
					s.addInPlace(m);
				
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
				
				final List<Smoother> ret = new ArrayList<>();
				for (int i = 1; i < fluid.refinements + 1; i++)
				{
					ret.add(new BSSmoother4(1, 1, this,
					                        spaces.get(i)
					                              .getVelocitySize()));
				}
				verbose = true;
				return ret;
			}
			
			@Override
			public void applyZeroBoundaryConditions(final TaylorHoodSpace space, final MutableVector vector)
			{
				final Int2DoubleMap nodeValues = fluid.getDirichletNodeValuesForSpace(space, t);
				nodeValues.forEach((node, val) -> vector.set(0, node));
			}
			
			@Override
			public AMGPreconditionerSpace<ContinuousTPFEVectorSpace, TPCell, TPFace, ContinuousTPVectorFunction,
				CoordinateVector, CoordinateMatrix, CoordinateTensor> getVelocityAMG()
			{
				final var me = this;
				return new AMGPreconditionerSpace<>(DLMLagrangeAMGSolver.this.fluid.refinements,
				                                    DLMLagrangeAMGSolver.this.fluid.polynomialDegree + 1)
				{
					@Override
					public List<ContinuousTPFEVectorSpace> createSpaces(final int refinements)
					{
						final List<ContinuousTPFEVectorSpace> ret = new ArrayList<>();
						for (int i = 0; i < refinements + 1; i++)
						{
							final ContinuousTPFEVectorSpace subSpace
								= new ContinuousTPFEVectorSpace(fluid.startCoordinates,
								                                fluid.endCoordinates,
								                                fluid.coarsestCells.mul(
									                                (int) Math.pow(
										                                2,
										                                i)));
							subSpace.assembleCells();
							subSpace.assembleFunctions(fluid.polynomialDegree + 1);
							ret.add(subSpace);
						}
						return ret;
					}
					
					@Override
					protected @NotNull SparseMatrix buildProlongationMatrix(
						final ContinuousTPFEVectorSpace coarse,
						final ContinuousTPFEVectorSpace fine)
					{
						return me.prolongationMatrices
							.get(me.getLevelFromVelocitySize(coarse.getShapeFunctions()
							                                       .size()))
							.slice(new IntCoordinates(0, 0),
							       new IntCoordinates(fine.getShapeFunctions()
							                              .size(),
							                          coarse.getShapeFunctions()
							                                .size()));
					}
					
					@Override
					public Tuple2<SparseMatrix, DenseVector> createSystem(final ContinuousTPFEVectorSpace continuousTPFEVectorSpace)
					{
						return new Tuple2<>(((SparseMatrix) me.getFinestSystem())
							                    .slice(new IntCoordinates(0, 0),
							                           new IntCoordinates(
								                           continuousTPFEVectorSpace.getShapeFunctions()
								                                                    .size(),
								                           continuousTPFEVectorSpace.getShapeFunctions()
								                                                    .size()
							                           )),
						                    new DenseVector(continuousTPFEVectorSpace.getShapeFunctions()
						                                                             .size()));
					}
					
					@Override
					public List<Smoother> createSmoothers()
					{
						final List<Smoother> ret = new ArrayList<>();
						for (int i = 1; i < fluid.refinements + 1; i++)
						{
							ret.add(new ForwardBackwardGaussSeidelSmoother(20,
							                                               systems.get(i)));
						}
						return ret;
					}
					
					@Override
					public void applyZeroBoundaryConditions(final ContinuousTPFEVectorSpace continuousTPFEVectorSpace,
					                                        final MutableVector vector)
					{
						final int level = me.getLevelFromVelocitySize(vector.getLength());
						
						final Int2DoubleMap nodeValues =
							fluid.getDirichletNodeValuesForSpace(me.spaces.get(level),
							                                     t);
						nodeValues.forEach((node, val) ->
						                   {
							                   if (node < vector.getLength())
								                   vector.set(0, node);
						                   });
					}
				};
			}
			
			private int getLevelFromVelocitySize(final int velocitySize)
			{
				for (int i = 0; i < spaces.size(); i++)
					if (spaces.get(i)
					          .getVelocitySize() == velocitySize)
						return i;
				throw new IllegalArgumentException("velocitySize not valid");
			}
		};
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
			schur = new RichardsonSchur(systemMatrix,
			                            null);
		else
			schur.resetOffDiagonals(systemMatrix);
		final AMGPreconditionerSpace<?, ?, ?, ?, ?, ?, ?> mg = create_space(particleStates, dt, t, fluidState);
		schur.preconditioner = mg;
		return schur.mvMul(rhs);
	}
	
	abstract class DLMAMG
		extends AMGPreconditionerSpace<TaylorHoodSpace, TPCell, TPFace, QkQkFunction, MixedValue, MixedGradient,
		MixedHessian>
	{
		
		public DLMAMG(final int refinements, final int polynomialDegree)
		{
			super(refinements, polynomialDegree);
		}
		
		abstract public AMGPreconditionerSpace<ContinuousTPFEVectorSpace, TPCell, TPFace, ContinuousTPVectorFunction,
			CoordinateVector, CoordinateMatrix, CoordinateTensor> getVelocityAMG();
	}
}
