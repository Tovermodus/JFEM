package dlm;

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

public class DLMLagrangeMGSolver
	extends DLMSolver
{
	public double smootherOmega = 1;
	public int smootherRuns = 7;
	
	PreconditionedIterativeImplicitSchur schur;
	MultiGridFluid fluid;
	final List<Particle> particles;
	
	public DLMLagrangeMGSolver(final MultiGridFluid f, final List<Particle> particles)
	{
		fluid = f;
		this.particles = particles;
	}
	
	@NotNull
	private MGPreconditionerSpace<TaylorHoodSpace, TPCell, TPFace, QkQkFunction, MixedValue, MixedGradient,
		MixedHessian> create_space(
		final List<ParticleIterate> particleStates, final double dt, final double t,
		final FluidIterate iterate)
	{
		return new MGPreconditionerSpace<>(fluid.refinements, fluid.polynomialDegree)
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
				velocity = new MixedTPFESpaceFunction<>(space.getShapeFunctionMap(),
				                                        restrictedIterate).getVelocityFunction();
				final FluidSystem fs = fluid.getFluidSystemForSpace(space, velocity, 0,
				                                                    restrictedIterate);
				final var blockRhs = Fluid.getBlockRhs(fs, dt);
				System.out.println("symmetry " + blockRhs._1.sub(blockRhs._1.transpose())
				                                            .absMaxElement());
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
				{
					s.addInPlace(m);
					System.out.println("symmetry " + s.sub(s.transpose())
					                                  .absMaxElement());
				}
				final Int2DoubleMap nodeValues = fluid.getDirichletNodeValuesForSpace(space, t);
				nodeValues.forEach((node, val) ->
				                   {
					                   s.deleteColumn(node);
					                   s.deleteRow(node);
					                   s.set(1, node, node);
				                   });
				System.out.println("symmetry " + s.sub(s.transpose())
				                                  .absMaxElement());
				return new Tuple2<>(s, new DenseVector(n));
			}
			
			@Override
			public List<Smoother> createSmoothers()
			{
				
				final List<Smoother> ret = new ArrayList<>();
				for (int i = 1; i < fluid.refinements + 1; i++)
				{
					ret.add(new BSSmoother2(smootherRuns, smootherOmega,
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
				nodeValues.forEach((node, val) ->
					                   vector.set(0, node));
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
			schur = new PreconditionedIterativeImplicitSchur(systemMatrix,
			                                                 null);
		else
			schur.resetOffDiagonals(systemMatrix);
		final MGPreconditionerSpace<?, ?, ?, ?, ?, ?, ?> mg = create_space(particleStates, dt, t, fluidState);
		schur.preconditioner = mg;
		return schur.mvMul(rhs);
	}
}
