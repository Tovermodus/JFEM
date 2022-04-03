import basic.*;
import distorted.*;
import distorted.geometry.DistortedCell;
import dlm.HyperbolicCartesianDistorted;
import io.vavr.Function2;
import io.vavr.Tuple2;
import linalg.*;
import mixed.*;
import multigrid.MGPReconditioner2;
import multigrid.MGSpace;
import multigrid.RichardsonSmoother;
import multigrid.Smoother;
import org.jetbrains.annotations.NotNull;
import tensorproduct.ContinuousTPShapeFunction;
import tensorproduct.ContinuousTPVectorFunction;
import tensorproduct.TPVectorCellIntegral;
import tensorproduct.TPVectorCellIntegralOnCell;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.concurrent.ConcurrentSkipListMap;
import java.util.function.Function;
import java.util.function.Predicate;
import java.util.function.ToDoubleFunction;

public class DLMRingPrecondMG
	extends HyperbolicCartesianDistorted<QkQkFunction, ContinuousTPShapeFunction, ContinuousTPVectorFunction>
{
	
	final Map<CoordinateVector, CoordinateVector> velocityValues;
	final Map<CoordinateVector, Double> pressureValues;
	final Map<CoordinateVector, CoordinateVector> XValues;
	final Map<CoordinateVector, CoordinateVector> XPrimeValues;
	final Map<CoordinateVector, CoordinateVector> UXValues;
	final List<CoordinateVector> plotPoints;
	VectorMultiplyable precond;
	private double lastPrecondTime;
	
	public DLMRingPrecondMG(final double dt,
	                        final int timeSteps,
	                        final TaylorHoodSpace backgroundSpace,
	                        final List<DistortedGridSpace<DistortedVectorShapeFunction, CoordinateVector, CoordinateMatrix, CoordinateTensor>> particleSpaces)
	{
		super(dt, timeSteps, backgroundSpace, particleSpaces);
		velocityValues = new ConcurrentSkipListMap<>();
		pressureValues = new ConcurrentSkipListMap<>();
		XValues = new ConcurrentSkipListMap<>();
		XPrimeValues = new ConcurrentSkipListMap<>();
		UXValues = new ConcurrentSkipListMap<>();
		plotPoints = backgroundSpace.generatePlotPoints(30);
	}
	
	@Override
	protected DenseVector getAdditionalRhs()
	{
		return new DenseVector(getSystemSize());
	}
	
	@Override
	protected Function2<Matrix, Vector, MutableVector> getSolver()
	{
		final IterativeSolver it = new IterativeSolver();
		it.gm.printProgress = true;
		it.showProgress = true;
		it.showInterrupt = false;
		return (A, b) ->
		{
			System.out.println(A.getClass());
			System.out.println(A.getRows());
			//final DenseMatrix Ad = new DenseMatrix(A);
			System.out.println(A.getRows());
			return new DenseVector(it.solvePGMRES(A, getPrecond(A), b, 1e-6));
		};
	}
	
	PreconditionedIterativeImplicitSchur ies;
	
	public VectorMultiplyable calculatePrecond(final Matrix A)
	{
		final int[] blocks = new int[1 + particleSpaces.size()];
		for (int i = 0; i < blocks.length; i++)
			blocks[i] = blockStarts[2 * i];
		System.out.println("copy");
		final BlockSparseMatrix schurShape = new BlockSparseMatrix(A, blocks);
		System.out.println("done");
		System.out.println("Generaate Schur Solver");
		if (ies == null)
		{
			PlotWindow.addPlot(new MatrixPlot(schurShape.getBlockMatrix(0, 0)));
			ies = new PreconditionedIterativeImplicitSchur(schurShape, getMGPreconditioner());
			System.out.println("ies built");
		} else
		{
			ies.getSchurBlock()
			   .overrideBy(schurShape.getBlockMatrix(0, 0));
			PlotWindow.addPlot(new MatrixPlot(schurShape, "block"));
			for (int i = 0; i < particleSpaces.size(); i++)
			{
				if (schurShape.getBlockMatrix(0, i + 1) != null)
					ies.getTopBlock(i)
					   .overrideBy(schurShape.getBlockMatrix(0, i + 1));
				if (schurShape.getBlockMatrix(i + 1, 0) != null)
					ies.getLeftBlock(i)
					   .overrideBy(schurShape.getBlockMatrix(i + 1, 0));
			}
		}
		return ies;
	}
	
	private VectorMultiplyable getMGPreconditioner()
	{
		final var mg = new MGSpace<TaylorHoodSpace, TPCell, TPFace, QkQkFunction, MixedValue, MixedGradient,
			MixedHessian>(backgroundRefines, backgroundDegree)
		{
			@Override
			public List<TaylorHoodSpace> createSpaces(final int refinements)
			{
				
				final ArrayList<TaylorHoodSpace> ret = new ArrayList<>();
				for (int i = 0; i < refinements + 1; i++)
				{
					ret.add(getTaylorHoodSpace(i));
				}
				return ret;
			}
			
			@Override
			public Tuple2<VectorMultiplyable, DenseVector> createSystem(final TaylorHoodSpace space)
			{
				final SparseMatrix s = new SparseMatrix(space.getShapeFunctionMap()
				                                             .size(),
				                                        space.getShapeFunctionMap()
				                                             .size());
				space.writeCellIntegralsToMatrix(getBackgroundMassIntegrals(), s);
				space.writeCellIntegralsToMatrix(getSemiImplicitBackgroundIntegrals(getUp().getVelocityFunction()),
				                                 s);
				s.mulInPlace(1. / dt);
				final DenseVector d = new DenseVector(space.getShapeFunctionMap()
				                                           .size());
				space.writeCellIntegralsToMatrix(getBackgroundIntegrals(), s);
				final VectorFunction velocityBoundaryFunction
					= VectorFunction.fromLambda(velocityBoundaryValues(),
					                            2, 2);
				space.writeBoundaryValuesTo(new ComposedMixedFunction(velocityBoundaryFunction),
				                            getDirichletBoundary(),
				                            (f, function) -> function.hasVelocityFunction(),
				                            s,
				                            d);
				final Optional<QkQkFunction> firstPressure =
					space.getShapeFunctionMap()
					     .values()
					     .stream()
					     .filter(ComposeMixedShapeFunction::hasPressureFunction)
					     .findAny();
				firstPressure.ifPresent(st -> space.overWriteValue(st.getGlobalIndex(),
				                                                   0,
				                                                   s,
				                                                   d));
				System.out.println(s.transpose()
				                    .sub(s)
				                    .absMaxElement());
				return new Tuple2<>(s, d);
			}
			
			@Override
			public List<Smoother> createSmoothers()
			{
				PlotWindow.addPlot(new MatrixPlot((SparseMatrix) finest_system, "precond"));
				final ArrayList<Smoother> ret = new ArrayList<>();
				for (int i = 1; i < spaces.size(); i++)
				{
					ret.add(new RichardsonSmoother(0.2, 20));//, d.getInvertedDiagonalMatrix()));
				}
				return ret;
			}
		};
		System.out.println("mg built");
		return new MGPReconditioner2(mg);
	}
	
	public VectorMultiplyable getPrecond(final Matrix A)
	{
		//if (precond == null || time > lastPrecondTime + dt * 20)
		{
			precond = calculatePrecond(A);
			System.out.println("calculated");
			lastPrecondTime = getTime();
		}
		return precond;
	}
	
	@Override
	protected void postIterationCallback()
	{
//		if (precond == null)
//		{
//			precond = getPrecond();
//		}
		final MixedFunctionOnCells<TPCell, TPFace> Up = getUp();
		velocityValues.putAll(Up.velocityValuesInPointsAtTime(plotPoints, getTime()));
		pressureValues.putAll(Up.pressureValuesInPointsAtTime(plotPoints, getTime()));
		final DistortedVectorFESpaceFunction X =
			new DistortedVectorFESpaceFunction(particleSpaces.get(0)
			                                                 .getShapeFunctionMap(),
			                                   getParticleIterate(getCurrentIterate(), 0));
		XValues.putAll(X.valuesInPointsAtTime(plotPoints, getTime()));
		final DistortedVectorFESpaceFunction XPrime =
			new DistortedVectorFESpaceFunction(particleSpaces.get(0)
			                                                 .getShapeFunctionMap(),
			                                   getParticleIterate(getCurrentIterate().sub(getLastIterate())
			                                                                         .mul(1. / dt), 0));
		XPrimeValues.putAll(XPrime.valuesInPointsAtTime(plotPoints, getTime()));
		final DistortedVectorFunctionOnCells UX =
			DistortedVectorFunctionOnCells.concatenate(Up.getVelocityFunction(), X);
		UXValues.putAll(UX.valuesInPointsAtTime(plotPoints, getTime()));
		System.out.println("ITeration at time " + getTime() + " out of " + dt * timeSteps + " is finished");
	}
	
	@Override
	protected List<CellIntegral<TPCell, QkQkFunction>> getBackgroundIntegrals()
	{
		final TPVectorCellIntegral<ContinuousTPVectorFunction> symGrad =
			new TPVectorCellIntegral<>(1, TPVectorCellIntegral.SYM_GRAD);
		final MixedCellIntegral<TPCell, ContinuousTPShapeFunction, ContinuousTPVectorFunction, QkQkFunction>
			symGradMixed = MixedTPCellIntegral.fromVelocityIntegral(symGrad);
		final MixedCellIntegral<TPCell, ContinuousTPShapeFunction, ContinuousTPVectorFunction, QkQkFunction>
			divValue
			= new MixedTPCellIntegral<>(ScalarFunction.constantFunction(-1), MixedTPCellIntegral.DIV_VALUE);
		return List.of(symGradMixed, divValue);
	}
	
	@Override
	protected List<CellIntegral<TPCell, QkQkFunction>> getSemiImplicitBackgroundIntegrals(final VectorFunctionOnCells<TPCell, TPFace> u)
	{
		
		final VectorFunctionOnCells<TPCell, TPFace> semiImplicitWeight1 =
			VectorFunctionOnCells.fromLambda((x) -> u.value(x)
			                                         .mul(1. / 2),
			                                 (x, cell) -> u.valueInCell(x, cell)
			                                               .mul(1. / 2), 2, 2);
		final VectorFunctionOnCells<TPCell, TPFace> semiImplicitWeight2 =
			VectorFunctionOnCells.fromLambda((x) -> u.value(x)
			                                         .mul(-1. / 2),
			                                 (x, cell) -> u.valueInCell(x, cell)
			                                               .mul(-1. / 2), 2, 2);
		final TPVectorCellIntegralOnCell<ContinuousTPVectorFunction> convection1 =
			new TPVectorCellIntegralOnCell<>(semiImplicitWeight1, TPVectorCellIntegral.GRAD_VALUE);
		final TPVectorCellIntegralOnCell<ContinuousTPVectorFunction> convection2 =
			new TPVectorCellIntegralOnCell<>(semiImplicitWeight2, TPVectorCellIntegral.VALUE_GRAD);
		final MixedCellIntegral<TPCell, ContinuousTPShapeFunction, ContinuousTPVectorFunction,
			QkQkFunction> mixedConvection1 = MixedCellIntegral.fromVelocityIntegral(convection1);
		final MixedCellIntegral<TPCell, ContinuousTPShapeFunction, ContinuousTPVectorFunction,
			QkQkFunction> mixedConvection2 = MixedCellIntegral.fromVelocityIntegral(convection2);
		return List.of(mixedConvection1, mixedConvection2);
	}
	
	@Override
	protected List<CellIntegral<DistortedCell, DistortedVectorShapeFunction>> getParticleIntegrals(final int particleId)
	{
		final double lameLambda;
		final double lameMu;
		if (particleId == 0)
		{
			lameLambda = 1;
			lameMu = 100;
		} else if (particleId == 1)
		{
			lameLambda = 1;
			lameMu = 100;
		} else// (particleId == 1)
		{
			lameLambda = 0.1;
			lameMu = 1;
		}
		final DistortedVectorCellIntegral lambdaIntegral = new DistortedVectorCellIntegral(lameLambda,
		                                                                                   DistortedVectorCellIntegral.TRACE_SYM);
		final DistortedVectorCellIntegral muIntegral = new DistortedVectorCellIntegral(lameMu * 2,
		                                                                               DistortedVectorCellIntegral.SYM_SYM);
		return List.of(lambdaIntegral, muIntegral);
	}
	
	@Override
	protected List<CellIntegral<DistortedCell, DistortedVectorShapeFunction>> getParticleLagrangeIntegrals(final int particleId)
	{
		final DistortedVectorCellIntegral h1Integral = new DistortedVectorCellIntegral(-1,
		                                                                               DistortedVectorCellIntegral.H1);
		return List.of(h1Integral);
	}
	
	@Override
	protected List<RightHandSideIntegral<DistortedCell, DistortedVectorShapeFunction>> getBackgroundLagrangeIntegrals(
		final QkQkFunction function,
		final DistortedVectorFunctionOnCells X,
		final int particleId)
	{
		final DistortedVectorDistortedRightHandSideIntegral lagrange =
			new DistortedVectorDistortedRightHandSideIntegral(DistortedVectorFunctionOnCells.concatenate(
				function.getVelocityFunction(), X),
			                                                  DistortedVectorDistortedRightHandSideIntegral.H1);
		return List.of(lagrange);
	}
	
	@Override
	protected List<CellIntegral<TPCell, QkQkFunction>> getBackgroundMassIntegrals()
	{
		final TPVectorCellIntegral<ContinuousTPVectorFunction> mass =
			new TPVectorCellIntegral<>(1, TPVectorCellIntegral.VALUE_VALUE);
		final MixedCellIntegral<TPCell, ContinuousTPShapeFunction, ContinuousTPVectorFunction, QkQkFunction>
			mixedMass
			= MixedCellIntegral.fromVelocityIntegral(mass);
		return List.of(mixedMass);
	}
	
	@Override
	protected List<CellIntegral<DistortedCell, DistortedVectorShapeFunction>> getParticleMassIntegrals(
		final int particleId)
	{
		final double density;
		if (particleId == 0)
			density = 10;
		else if (particleId == 1)
			density = 10;
		else
			density = 1;
		final DistortedVectorCellIntegral massIntegral = new DistortedVectorCellIntegral(density,
		                                                                                 DistortedVectorCellIntegral.VALUE_VALUE);
		return List.of(massIntegral);
	}
	
	@Override
	protected Function<CoordinateVector, CoordinateVector> getInitialFluidVelocity()
	{
		return x -> new CoordinateVector(x.getLength());
	}
	
	@Override
	protected ToDoubleFunction<CoordinateVector> getInitialPressure()
	{
		return x -> 0;
	}
	
	@Override
	protected Function<CoordinateVector, CoordinateVector> getInitialParticleVelocity(final int particleId)
	{
		System.out.println("initial particle velo" + particleId);
		if (particleId == 0)
			return x -> CoordinateVector.fromValues(-1.0, .3);
		if (particleId == 1)
			return x -> CoordinateVector.fromValues(1.0, .3);
		if (particleId == 2)
			return x -> CoordinateVector.fromValues(0, -.5);
		return x -> new CoordinateVector(x.getLength());
	}
	
	@Override
	protected Function<CoordinateVector, CoordinateVector> velocityBoundaryValues()
	{
		return x -> CoordinateVector.getUnitVector(2, 0)
		                            .mul(0.1);
	}
	
	@Override
	protected Predicate<TPFace> getDirichletBoundary()
	{
		return f -> f.center()
		             .y() == 0 || f.center()
		                           .y() == 1;
	}
	
	private static final int backgroundRefines = 3;
	private static final int backgroundDegree = 1;
	
	public static void main(final String[] args)
	{
		final TaylorHoodSpace backGround = getTaylorHoodSpace(backgroundRefines);
		final DistortedGridSpace<DistortedVectorShapeFunction, CoordinateVector, CoordinateMatrix, CoordinateTensor>
			particle1 =
			new CircleVectorSpace(CoordinateVector.fromValues(0.6, 0.5),
			                      0.05,
			                      1);
		final DistortedGridSpace<DistortedVectorShapeFunction, CoordinateVector, CoordinateMatrix, CoordinateTensor>
			particle0 =
			new CircleVectorSpace(CoordinateVector.fromValues(0.4, 0.5),
			                      0.05, 1);
		final DistortedGridSpace<DistortedVectorShapeFunction, CoordinateVector, CoordinateMatrix, CoordinateTensor>
			particle2 =
			new RingVectorSpace(CoordinateVector.fromValues(0.5, 0.5),
			                    0.25,
			                    0.255,
			                    0);
		backGround.assembleCells();
		backGround.assembleFunctions(backgroundDegree);
		particle0.assembleCells();
		particle0.assembleFunctions(1);
		particle1.assembleCells();
		particle1.assembleFunctions(1);
		particle2.assembleCells();
		particle2.assembleFunctions(1);
		final DLMRingPrecondMG dlmElast2 = new DLMRingPrecondMG(0.02,
		                                                        10,
		                                                        backGround,
		                                                        List.of(particle0, particle1, particle2));
		dlmElast2.loop();
		final MixedPlot2DTime UpPlot0 = new MixedPlot2DTime(dlmElast2.pressureValues,
		                                                    dlmElast2.velocityValues,
		                                                    30,
		                                                    "Up1");
		final Overlay ball1 = new DistortedOverlay(dlmElast2.plotPoints,
		                                           particle0,
		                                           dlmElast2.getIterateHistory()
		                                                    .slice(new IntCoordinates(0,
		                                                                              backGround.getShapeFunctionMap()
		                                                                                        .size()),
		                                                           new IntCoordinates(dlmElast2.getIterateHistory()
		                                                                                       .getRows(),
		                                                                              backGround.getShapeFunctionMap()
		                                                                                        .size()
			                                                                              + particle0.getShapeFunctionMap()
			                                                                                         .size())),
		                                           5);
		UpPlot0.addOverlay(ball1);
		
		PlotWindow.addPlot(UpPlot0);
		final MixedPlot2DTime UpPlot = new MixedPlot2DTime(dlmElast2.pressureValues,
		                                                   dlmElast2.velocityValues,
		                                                   30,
		                                                   "Up1");
		final Overlay ball2 = new DistortedOverlay(dlmElast2.plotPoints,
		                                           particle1,
		                                           dlmElast2.getIterateHistory()
		                                                    .slice(new IntCoordinates(0,
		                                                                              backGround.getShapeFunctionMap()
		                                                                                        .size()
			                                                                              + 2 * particle0.getShapeFunctionMap()
			                                                                                             .size()),
		                                                           new IntCoordinates(dlmElast2.getIterateHistory()
		                                                                                       .getRows(),
		                                                                              backGround.getShapeFunctionMap()
		                                                                                        .size()
			                                                                              + 2 * particle0.getShapeFunctionMap()
			                                                                                             .size()
			                                                                              + particle1.getShapeFunctionMap()
			                                                                                         .size())),
		                                           5);
		UpPlot.addOverlay(ball2);
		PlotWindow.addPlot(UpPlot);
		final MixedPlot2DTime UpPlot1 = new MixedPlot2DTime(dlmElast2.pressureValues,
		                                                    dlmElast2.velocityValues,
		                                                    30,
		                                                    "Up2");
		final Overlay ring = new DistortedOverlay(dlmElast2.plotPoints,
		                                          particle2,
		                                          dlmElast2.getIterateHistory()
		                                                   .slice(new IntCoordinates(0,
		                                                                             backGround.getShapeFunctionMap()
		                                                                                       .size()
			                                                                             + 2 * particle0.getShapeFunctionMap()
			                                                                                            .size()
			                                                                             + 2 * particle1.getShapeFunctionMap()
			                                                                                            .size()),
		                                                          new IntCoordinates(dlmElast2.getIterateHistory()
		                                                                                      .getRows(),
		                                                                             backGround.getShapeFunctionMap()
		                                                                                       .size()
			                                                                             + 2 * particle0.getShapeFunctionMap()
			                                                                                            .size()
			                                                                             + 2 * particle1.getShapeFunctionMap()
			                                                                                            .size()
			                                                                             + particle2.getShapeFunctionMap()
			                                                                                        .size())),
		                                          5);
		UpPlot1.addOverlay(ring);
		UpPlot0.addOverlay(ring);
		UpPlot0.addOverlay(ball2);
		final Overlay coordinateSys = new CoordinateSystemOverlay(CoordinateVector.fromValues(0, 0),
		                                                          CoordinateVector.fromValues(1, 1));
		UpPlot0.addOverlay(coordinateSys);
		PlotWindow.addPlot(UpPlot1);
	}
	
	@NotNull
	private static TaylorHoodSpace getTaylorHoodSpace(final int refines)
	{
		final int mul = (int) Math.pow(2, refines);
		return new TaylorHoodSpace(CoordinateVector.fromValues(0, 0),
		                           CoordinateVector.fromValues(1, 1),
		                           new IntCoordinates(2 * mul, 2 * mul));
	}
}
