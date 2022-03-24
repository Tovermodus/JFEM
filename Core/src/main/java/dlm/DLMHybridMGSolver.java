package dlm;

import basic.PlotWindow;
import basic.VectorFunctionOnCells;
import io.vavr.Tuple2;
import it.unimi.dsi.fastutil.ints.Int2DoubleMap;
import linalg.*;
import mixed.*;
import multigrid.MGPreconditionerSpace;
import multigrid.Smoother;
import org.jetbrains.annotations.NotNull;
import schwarz.DirectSolver;
import schwarz.MultiplicativeSubspaceCorrection;
import schwarz.SchwarzSmoother;
import schwarz.VankaSchwarz;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class DLMHybridMGSolver
	extends DLMSolver
{
	public final int smootherRuns;
	final int overlap;
	MultigridFixedpointSchur schur;
	MultiGridFluid fluid;
	final List<Particle> particles;
	final double omega;
	
	public DLMHybridMGSolver(final int smootherRuns,
	                         final int overlap,
	                         final MultiGridFluid f,
	                         final List<Particle> particles, final double omega)
	{
		this.smootherRuns = smootherRuns;
		this.overlap = overlap;
		fluid = f;
		this.particles = particles;
		this.omega = omega;
	}
	
	private void addSchurAMG(final MGPreconditionerSpace<TaylorHoodSpace, TPCell, TPFace, QkQkFunction, MixedValue, MixedGradient,
		MixedHessian> mg,
	                         final List<ParticleIterate> particleStates,
	                         final List<ParticleSystem> particleSystems)
	{
		SparseMatrix finerSchurContribution = new SparseMatrix(mg.finest_system.getVectorSize(),
		                                                       mg.finest_system.getTVectorSize());
		final List<Matrix> schurContributions =
			IntStream.range(0, particles.size())
			         .parallel()
			         .mapToObj(i -> new Tuple2<>(i, particleSystems.get(i).lagrangeBackgroundMatrix))
			         .map(e ->
			              {
				              final SparseMatrix block =
					              new SparseMatrix(mg.finest_system.getVectorSize(),
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
			finerSchurContribution.addInPlace(m);
		}
		((SparseMatrix) mg.getFinestSystem()).addInPlace(finerSchurContribution);
		for (int coarserSystemIndex = mg.systems.size() - 2; coarserSystemIndex >= 0; coarserSystemIndex--)
		{
			final SparseMatrix coarserSystemMatrix = (SparseMatrix) mg.getSystem(coarserSystemIndex);
			final SparseMatrix prolongationFromCoarser = mg.getProlongationMatrix(coarserSystemIndex);
			final Matrix coarserSchurContribution =
				prolongationFromCoarser.tmMul(finerSchurContribution.mmMul(prolongationFromCoarser));
			coarserSystemMatrix.addInPlace(coarserSchurContribution);
			finerSchurContribution = new SparseMatrix(coarserSchurContribution);
		}
	}
	
	private void fixBoundaryNodes(final MGPreconditionerSpace<TaylorHoodSpace, TPCell, TPFace, QkQkFunction, MixedValue, MixedGradient,
		MixedHessian> mg, final double t)
	{
		for (int i = 0; i < mg.systems.size(); i++)
		{
			final var space = mg.spaces.get(i);
			final SparseMatrix s = (SparseMatrix) mg.systems.get(i);
			final Int2DoubleMap nodeValues = fluid.getDirichletNodeValuesForSpace(space, t);
			nodeValues.forEach((node, val) ->
			                   {
				                   s.deleteColumn(node);
				                   s.deleteRow(node);
				                   s.set(1, node, node);
			                   });
		}
	}
	
	@NotNull
	private MGPreconditionerSpace<TaylorHoodSpace, TPCell, TPFace, QkQkFunction, MixedValue, MixedGradient,
		MixedHessian> create_space(
		final List<ParticleIterate> particleStates, final double dt, final double t,
		final FluidIterate iterate, final List<ParticleSystem> particleSystems)
	{
		return new MGPreconditionerSpace<>(fluid.refinements, fluid.polynomialDegree)
		{
			@Override
			public void postmoothcallback(final int level, final Vector guess, final Vector rhs)
			{
				if (level == 0)
				{
					final var function =
						new MixedTPFESpaceFunction<>(getSpace(level).getShapeFunctionMap(),
						                             guess);
					PlotWindow.addPlot(new MixedPlot2D(function,
					                                   getSpace(level).generatePlotPoints(40), 40
						, "coarse correction"));
				}
			}
			
			@Override
			public void precorrectioncallback(final int level, final Vector guess, final Vector rhs)
			{
				if (level == 1)
				{
					final Vector realCorrect =
						((SparseMatrix) getSystem(level)).solve(
							rhs.sub(getSystem(level).mvMul(guess)));
					final var function =
						new MixedTPFESpaceFunction<>(getSpace(level).getShapeFunctionMap(),
						                             realCorrect);
					PlotWindow.addPlot(new MixedPlot2D(function,
					                                   getSpace(level).generatePlotPoints(40), 40
						, "real pre correction"));
				}
			}
			
			@Override
			public void correctioncallback(final int level, final Vector correction, final Vector rhs)
			{
				if (level == 1)
				{
					final var function =
						new MixedTPFESpaceFunction<>(getSpace(level).getShapeFunctionMap(),
						                             correction);
					PlotWindow.addPlot(new MixedPlot2D(function,
					                                   getSpace(level).generatePlotPoints(40), 40
						, "actual correction"));
				}
			}
			
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
				Vector current = iterate.current;
				if (current == null)
					current = new DenseVector(getFinestSpace().getShapeFunctions()
					                                          .size());
				velocity = new MixedTPFESpaceFunction<>(getFinestSpace().getShapeFunctionMap(),
				                                        current).getVelocityFunction();
				
				final FluidSystem fs;
				if (space == getFinestSpace())
					fs = fluid.finestFliudSystem;
				else
					fs = fluid.getFluidSystemForSpace(space, velocity, 0,
					                                  new DenseVector(n));
				final var blockRhs = Fluid.assembleBlockRhs(fs, dt);
				final SparseMatrix s = new SparseMatrix(blockRhs._1);
				return new Tuple2<>(s, new DenseVector(n));
			}
			
			@Override
			public List<Smoother> createSmoothers()
			{
				
				addSchurAMG(this, particleStates, particleSystems);
				fixBoundaryNodes(this, t);
				final List<Smoother> ret = new ArrayList<>();
				for (int i = 1; i < fluid.refinements + 1; i++)
				{
					final IntCoordinates partitions =
						new IntCoordinates((int) (fluid.coarsestCells.get(0) * Math.pow(2,
						                                                                i - 2)),
						                   (int) (fluid.coarsestCells.get(1) * Math.pow(2,
						                                                                i - 2)));
					System.out.println(partitions);

//					final ColoredCartesianSchwarz<QkQkFunction> schwarz
//						= new ColoredCartesianSchwarz<>((SparseMatrix) getSystem(i),
//						                                getSpace(i),
//						                                partitions, overlap,
//						                                new DirectSolver(), omega);
//					final CartesianUpFrontSchwarz<QkQkFunction> schwarz
//						= new CartesianUpFrontSchwarz<>((SparseMatrix) getSystem(i),
//						                                getSpace(i),
//						                                partitions, overlap,
//						                                new AdditiveSubspaceCorrection<>(0.1),
//						                                new DirectSolver());
					final VankaSchwarz schwarz = new VankaSchwarz((SparseMatrix) getSystem(i),
					                                              getSpace(i),
					                                              new MultiplicativeSubspaceCorrection<>(),
					                                              new DirectSolver());
					ret.add(new SchwarzSmoother(smootherRuns, schwarz));
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
			schur = new MultigridFixedpointSchur(systemMatrix, null);//PreconditionedIterativeImplicitSchur
			// (systemMatrix, null);
		else
			schur.resetOffDiagonals(systemMatrix);
		final MGPreconditionerSpace<TaylorHoodSpace, TPCell, TPFace, QkQkFunction, MixedValue, MixedGradient,
			MixedHessian> mg = create_space(particleStates, dt, t, fluidState, particleSystems);
		mg.cycles = 2;
		schur.mg = mg;
//		final DirectSchur ds = new DirectSchur(systemMatrix);
//		System.out.println("Schur max eigenvalue " + ds.getSchurComplement()
//		                                               .powerIterationNonSymmetric());
//
//		System.out.println("Fluid max eigenvalue " + ds.getSchurBlock()
//		                                               .powerIterationNonSymmetric());
//
//		System.out.println("Schur min eigenvalue " + 1. / ds.getSchurComplement()
//		                                                    .inverse()
//		                                                    .powerIterationNonSymmetric());
//
//		System.out.println("Fluid min eigenvalue " + 1. / schur.getSchurBlock()
//		                                                       .inverse()
//		                                                       .powerIterationNonSymmetric());
		return schur.mvMul(rhs);
	}
}
