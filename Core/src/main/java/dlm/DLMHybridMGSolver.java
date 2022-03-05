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

public class DLMHybridMGSolver
	extends DLMSolver
{
	public double smootherOmega = 1;
	public int smootherRuns = 1;
	
	RichardsonSchur schur;
	MultiGridFluid fluid;
	final List<Particle> particles;
	
	public DLMHybridMGSolver(final MultiGridFluid f, final List<Particle> particles)
	{
		fluid = f;
		this.particles = particles;
	}
	
	private void addSchurAMG(final MGPreconditionerSpace<TaylorHoodSpace, TPCell, TPFace, QkQkFunction, MixedValue, MixedGradient,
		MixedHessian> mg, final List<ParticleIterate> particleStates)
	{
		SparseMatrix finerSchurContribution = new SparseMatrix(mg.finest_system.getVectorSize(),
		                                                       mg.finest_system.getTVectorSize());
		final List<Matrix> schurContributions =
			IntStream.range(0, particles.size())
			         .parallel()
			         .mapToObj(i -> new Tuple2<>(i,
			                                     particles.get(i)
			                                              .buildLagrangeBackgroundMatrix(
				                                              mg.getFinestSpace(),
				                                              particleStates.get(i))))
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
		final FluidIterate iterate)
	{
		return new MGPreconditionerSpace<>(fluid.refinements, fluid.polynomialDegree)
		{
//			@Override
//			public void presmoothcallback(final int level, final Vector guess, final Vector rhs)
//			{
//				if (level == 1)
//				{
//					final Vector sol = ((SparseMatrix) getFinestSystem()).solve(rhs);
//					final MixedTPFESpaceFunction<QkQkFunction> solutionFunction =
//						new MixedTPFESpaceFunction<>(getFinestSpace().getShapeFunctionMap(),
//						                             sol);
//					PlotWindow.addPlot(new MixedPlot2D(solutionFunction,
//					                                   getFinestSpace().generatePlotPoints(50),
//					                                   50, "solution"));
//					final MixedTPFESpaceFunction<QkQkFunction> guessFUnction =
//						new MixedTPFESpaceFunction<>(getFinestSpace().getShapeFunctionMap(),
//						                             guess);
//					PlotWindow.addPlot(new MixedPlot2D(guessFUnction,
//					                                   getFinestSpace().generatePlotPoints(50),
//					                                   50, "presmoot"));
//				}
//			}
//
//			@Override
//			public void postmoothcallback(final int level, final Vector guess, final Vector rhs)
//			{
//				final MixedTPFESpaceFunction<QkQkFunction> guessFUnction =
//					new MixedTPFESpaceFunction<>(getFinestSpace().getShapeFunctionMap(),
//					                             guess);
//				PlotWindow.addPlot(new MixedPlot2D(guessFUnction,
//				                                   getFinestSpace().generatePlotPoints(50),
//				                                   50, "postsmooth"));
//			}
//
//			@Override
//			public void precorrectioncallback(final int level, final Vector guess, final Vector rhs)
//			{
//				final MixedTPFESpaceFunction<QkQkFunction> guessFUnction =
//					new MixedTPFESpaceFunction<>(getFinestSpace().getShapeFunctionMap(),
//					                             guess);
//				PlotWindow.addPlot(new MixedPlot2D(guessFUnction,
//				                                   getFinestSpace().generatePlotPoints(50),
//				                                   50, "precorrect"));
//			}
//
//			@Override
//			public void postcorrectioncallback(final int level, final Vector guess, final Vector rhs)
//			{
//				final MixedTPFESpaceFunction<QkQkFunction> guessFUnction =
//					new MixedTPFESpaceFunction<>(getFinestSpace().getShapeFunctionMap(),
//					                             guess);
//				PlotWindow.addPlot(new MixedPlot2D(guessFUnction,
//				                                   getFinestSpace().generatePlotPoints(50),
//				                                   50, "postcorrect"));
//			}
			
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
				return new Tuple2<>(s, new DenseVector(n));
			}
			
			@Override
			public List<Smoother> createSmoothers()
			{
				
				addSchurAMG(this, particleStates);
				fixBoundaryNodes(this, t);
				final List<Smoother> ret = new ArrayList<>();
				for (int i = 1; i < fluid.refinements + 1; i++)
				{
//					final int vel_size = getSpace(i).getVelocitySize();
//					final int tot_size = getSpace(i).getShapeFunctions()
//					                                .size();
//					final SparseMatrix[] blocks =
//						((SparseMatrix) getSystem(i)).partition(new IntCoordinates(
//							vel_size,
//							vel_size));
//					final SparseMatrix B = blocks[2];
//					final SparseMatrix A = blocks[0];
//					final Map<IntCoordinates, SparseMatrix> blockMap = new HashMap<>();
//					blockMap.put(new IntCoordinates(0, 0), //A);
//					             SparseMatrix.identity(vel_size)
//					                         .mul(1000));
////                                                     A.getDiagonalMatrix()
////                                                      .mul(200));
//					blockMap.put(new IntCoordinates(0, vel_size),
//					             B.transpose());
//					blockMap.put(new IntCoordinates(vel_size, 0),
//					             B);
//					final BlockSparseMatrix p = new BlockSparseMatrix(blockMap,
//					                                                  tot_size,
//					                                                  tot_size);
//					final SparseMatrix prec = p.toSparse();
//					final Int2DoubleMap nodeValues =
//						fluid.getDirichletNodeValuesForSpace(getSpace(i), t);
//					nodeValues.forEach((node, val) ->
//					                   {
//						                   prec.deleteColumn(node);
//						                   prec.deleteRow(node);
//						                   prec.set(1, node, node);
//					                   });
//					final int firstPressure =
//						getSpace(i).getVelocitySize();
//					getSpace(i).overWriteValue(firstPressure, 0, prec,
//					                           new DenseVector(prec.getRows()));
//					final DenseMatrix pin = prec.inverse();
//					final SparseMatrix alphaIB = prec.sub((SparseMatrix) getSystem(i));
//					final Matrix pinAlphaIB = pin.mmMul(alphaIB);
//					PlotWindow.addPlot(new MatrixPlot(pinAlphaIB));
//					PlotWindow.addPlot(new MatrixPlot(pinAlphaIB.slice(new IntCoordinates(
//						                                                   0,
//						                                                   0)
//						, new IntCoordinates(vel_size, vel_size))));
////					System.out.println("MAX EIG errit" + pinAlphaIB.slice(new IntCoordinates(
////						                                                      0,
////						                                                      0)
////						                                               , new IntCoordinates(vel_size, vel_size))
////					                                               .powerIterationNonSymmetric());
//					ret.add(new RichardsonSmoother(0.1,
//					                               10,
//					                               pin));
					ret.add(new SIMPLEAMG(3, 0.02, getSpace(i).getVelocitySize(), this));
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
			schur = new RichardsonSchur(systemMatrix,
			                            null);
		else
			schur.resetOffDiagonals(systemMatrix);
		final MGPreconditionerSpace<TaylorHoodSpace, TPCell, TPFace, QkQkFunction, MixedValue, MixedGradient,
			MixedHessian> mg = create_space(particleStates, dt, t, fluidState);
		schur.preconditioner = mg;
		return schur.mvMul(rhs);
	}
}
