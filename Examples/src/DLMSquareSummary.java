import basic.*;
import linalg.*;
import mixed.*;
import tensorproduct.*;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

public class DLMSquareSummary
{
	final double dt;
	final int timeSteps;
	final double rhoF;
	final double rhoS;
	final double nu;
	final double kappa;
	TaylorHoodSpace eulerian;
	ContinuousTPFEVectorSpace lagrangian;
	SparseMatrix lagrangianBetaMass;
	SparseMatrix eulerianAlphaMass;
	SparseMatrix lagrangianDual;
	int nEulerian;
	int nLagrangian;
	int nTransfer;
	int eulerianPointsPerDimension = 30;
	int nEulerCells = 6;
	int nLagrangianCells = 5;
	List<CoordinateVector> eulerianPoints;
	private final int lagrangeDegree = 1;
	private final int eulerDegree = 1;
	private final String elastMethod = TPVectorCellIntegral.SYM_GRAD;
	
	public DLMSquareSummary()
	{
//		PerformanceArguments.PerformanceArgumentBuilder builder =
//			new PerformanceArguments.PerformanceArgumentBuilder();
//		builder.executeChecks = false;
//		builder.build();
		rhoF = 1;
		rhoS = 2;
		nu = 1;
		kappa = 100.000000;
		dt = 0.05;
		timeSteps = 20;
		initializeEulerian();
		initializeLagrangian();
		
		final SparseMatrix constantSystemMatrix = new SparseMatrix(nEulerian + nLagrangian + nTransfer,
		                                                           nEulerian + nLagrangian + nTransfer);
		SparseMatrix systemMatrix;
		DenseVector rightHandSide;
		final DenseVector constantRightHandSide = new DenseVector(nEulerian + nLagrangian + nTransfer);
		writeABf(constantSystemMatrix);
		System.out.println("abf");
		writeAs(constantSystemMatrix);
		System.out.println("as");
		
		writeCs(constantSystemMatrix);
		System.out.println("Cs");
		DenseVector currentIterate = new DenseVector(nEulerian + nLagrangian + nTransfer);
		DenseVector lastIterate = new DenseVector(nEulerian + nLagrangian + nTransfer);
		generateFirstEulerianIterates(currentIterate);
		generateFirstLagrangianIterates(currentIterate, lastIterate);
		System.out.println("firstIts");
		writeBoundaryValues(constantSystemMatrix, currentIterate, 0);
		writeBoundaryValues(constantSystemMatrix, lastIterate, -dt);
		System.out.println("bdrs");
		final Map<CoordinateVector, Double> pvals = getUp(currentIterate)
			.pressureValuesInPointsAtTime(eulerianPoints, -dt);
		final Map<CoordinateVector, CoordinateVector> vvals = getUp(currentIterate)
			.velocityValuesInPointsAtTime(eulerianPoints, -dt);
		VectorFunction X = getX(lastIterate);
		//final ScalarFunction circleIndicator = lagrangian.getIndicatorFunction();
		final VectorFunction finalX = X;
		final Map<CoordinateVector, Double> dvals = eulerianPoints
			.stream()
			.map(p ->
			     {
				     double val = 1;
				     CoordinateVector Xp = finalX.value(p);
				     if (Xp.euclidianNorm() == 0)
				     {
					     Xp = p;
					     val = 0;
				     }
				     return Map.entry(Xp.addCoordinate(-dt), val);
			     })
			.collect(Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue));
		final Map<CoordinateVector, CoordinateVector> derivVals = getX(lastIterate.mul(0))
			.valuesInPointsAtTime(eulerianPoints, -dt);
		final Map<CoordinateVector, CoordinateVector> uXvals = concatenateVelocityWithX(getUp(lastIterate),
		                                                                                lastIterate)
			.valuesInPointsAtTime(eulerianPoints, -dt);
		
		final DenseMatrix iterateHistory = new DenseMatrix(timeSteps + 1, nEulerian + nLagrangian + nTransfer);
		final DenseMatrix rhsHistory = new DenseMatrix(timeSteps + 1, nEulerian + nLagrangian + nTransfer);
		iterateHistory.addRow(lastIterate, 0);
		iterateHistory.addRow(currentIterate, 1);
		for (int i = 1; i < timeSteps; i++)
		{
			rightHandSide = new DenseVector(constantRightHandSide);
			systemMatrix = new SparseMatrix(constantSystemMatrix);
			pvals.putAll(getUp(currentIterate).pressureValuesInPointsAtTime(eulerianPoints, (i - 1) * dt));
			vvals.putAll(getUp(currentIterate).velocityValuesInPointsAtTime(eulerianPoints, (i - 1) * dt));
			X = getX(currentIterate);
			final VectorFunction finalX1 = X;
			final int finalI = i;
			dvals.putAll(eulerianPoints
				             .stream()
				             .map(p ->
				                  {
					                  double val = 1;
					                  CoordinateVector Xp = finalX1.value(p);
					                  if (Xp.euclidianNorm() == 0)
					                  {
						                  Xp = p;
						                  val = 0;
					                  }
					                  return Map.entry(Xp.addCoordinate((finalI - 1) * dt), val);
				                  })
				             .collect(Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue)));
			derivVals.putAll(getX(currentIterate.sub(lastIterate)
			                                    .mul(1. / dt))
				                 .valuesInPointsAtTime(eulerianPoints, (i - 1) * dt));
			uXvals.putAll(concatenateVelocityWithX(getUp(currentIterate), currentIterate)
				              .valuesInPointsAtTime(eulerianPoints, (i - 1) * dt));
			System.out.println("saved Iterate");
			writeF(rightHandSide, currentIterate);
			System.out.println("f");
			writeG(rightHandSide, currentIterate, lastIterate);
			System.out.println("g");
			writeD(rightHandSide, currentIterate);
			System.out.println("d");
			writeCf(systemMatrix, currentIterate);
			System.out.println("cf");
			writeBoundaryValues(systemMatrix, rightHandSide, i * dt);
			System.out.println("bdr");
			System.out.println(i + "th iteration");
			lastIterate = new DenseVector(currentIterate);
			rhsHistory.addRow(rightHandSide, i + 1);
			
			currentIterate = systemMatrix.solve(rightHandSide);
			System.out.println("newit" + getLagrangianIterate(currentIterate));
			iterateHistory.addRow(currentIterate, i + 1);
			System.out.println("solved");
		}
		PlotWindow.addPlot(new MatrixPlot(iterateHistory, "iterateHistory"));
		PlotWindow.addPlot(new MatrixPlot(rhsHistory, "rhsHistory"));
		final MixedPlot2DTime plot = new MixedPlot2DTime(pvals, vvals, eulerianPointsPerDimension, "VVals");
		PlotWindow.addPlot(new MixedPlot2DTime(pvals, uXvals, eulerianPointsPerDimension, "uXVals"));
		try
		{
			Thread.sleep(1000);
		} catch (final InterruptedException e)
		{
			e.printStackTrace();
		}
		PlotWindow.addPlot(new MixedPlot2DTime(pvals, derivVals, eulerianPointsPerDimension, "derivVals"));
		
		PlotWindow.addPlot(plot.addOverlay(new CartesianOverlay<>(eulerianPoints, lagrangian,
		                                                          iterateHistory.slice(new IntCoordinates(0,
		                                                                                                  nEulerian),
		                                                                               new IntCoordinates(
			                                                                               timeSteps + 1,
			                                                                               nEulerian + nLagrangian)),
		                                                          10)));
		final ScalarPlot2DTime disPLacementplot = new ScalarPlot2DTime(dvals, eulerianPointsPerDimension,
		                                                               "dVals");
		PlotWindow.addPlot(disPLacementplot);
	}
	
	public static void main(final String[] args)
	{
		new DLMSquareSummary();
	}
	
	public static VectorFunction initialVelocityValues()
	{
		return new VectorFunction()
		{
			@Override
			public int getRangeDimension()
			{
				return 2;
			}
			
			@Override
			public int getDomainDimension()
			{
				return 2;
			}
			
			@Override
			public CoordinateVector value(final CoordinateVector pos)
			{
				return CoordinateVector.fromValues(0, 0.0);
			}
		};
	}
	
	public static VectorFunction boundaryVelocityValues(final double t)
	{
		return new VectorFunction()
		{
			@Override
			public int getRangeDimension()
			{
				return 2;
			}
			
			@Override
			public int getDomainDimension()
			{
				return 3;
			}
			
			@Override
			public CoordinateVector value(final CoordinateVector pos)
			{
				return CoordinateVector.fromValues(Math.max(0, 0.01 * Math.min(0.3, t * 10)), 0);
			}
		};
	}
	
	public static VectorFunction initialDisplacementValues()
	{
		return new VectorFunction()
		{
			@Override
			public int getRangeDimension()
			{
				return 2;
			}
			
			@Override
			public int getDomainDimension()
			{
				return 2;
			}
			
			@Override
			public CoordinateVector value(final CoordinateVector pos)
			{
				return pos;
			}
		};
	}
	
	public static VectorFunction initialDisplacementDerivatives()
	{
		
		return new VectorFunction()
		{
			@Override
			public int getRangeDimension()
			{
				return 2;
			}
			
			@Override
			public int getDomainDimension()
			{
				return 2;
			}
			
			@Override
			public CoordinateVector value(final CoordinateVector pos)
			{
				return CoordinateVector.fromValues(-0.000, 0.0000);
			}
		};
	}
	
	public static ScalarFunction initialPressureValues()
	{
		return new ScalarFunction()
		{
			
			@Override
			public int getDomainDimension()
			{
				return 2;
			}
			
			@Override
			public Double value(final CoordinateVector pos)
			{
				return 0.;
			}
		};
	}
	
	public void initializeEulerian()
	{
		final CoordinateVector startCoordinates = CoordinateVector.fromValues(0, 0);
		final CoordinateVector endCoordinates = CoordinateVector.fromValues(1, 1);
		final IntCoordinates cellCounts = new IntCoordinates(nEulerCells, nEulerCells);
		eulerian = new TaylorHoodSpace(startCoordinates, endCoordinates, cellCounts);
		eulerian.assembleCells();
		eulerian.assembleFunctions(eulerDegree);
		nEulerian = eulerian.getShapeFunctionMap()
		                    .size();
		eulerianPoints = eulerian.generatePlotPoints(eulerianPointsPerDimension);
	}
	
	public void initializeLagrangian()
	{
		final CoordinateVector center = CoordinateVector.fromValues(0.5, 0.5);
		final double radius = 0.2;
		lagrangian = new ContinuousTPFEVectorSpace(center.sub(CoordinateVector.fromValues(0.2, 0.2)),
		                                           center.add(CoordinateVector.fromValues(0.2, 0.2)),
		                                           new IntCoordinates(nLagrangianCells, nLagrangianCells));
		lagrangian.assembleCells();
		lagrangian.assembleFunctions(lagrangeDegree);
		nLagrangian = lagrangian.getShapeFunctionMap()
		                        .size();
		nTransfer = nLagrangian;
	}
	
	public void generateFirstEulerianIterates(final DenseVector initialVector)
	{
		final MixedFunction initialVelocityFunction = new ComposedMixedFunction(initialVelocityValues());
		final MixedFunction initialPressureFunction = new ComposedMixedFunction(initialPressureValues());
		eulerian
			.getShapeFunctionMap()
			.forEach((key, value) -> initialVector.add(
				value.getNodeFunctional()
				     .evaluate(initialVelocityFunction), key));
		eulerian
			.getShapeFunctionMap()
			.forEach((key, value) -> initialVector.add(
				value.getNodeFunctional()
				     .evaluate(initialPressureFunction), key));
	}
	
	public void generateFirstLagrangianIterates(final DenseVector initialVector, final DenseVector previousVector)
	{
		lagrangian
			.getShapeFunctionMap()
			.forEach((key, value) -> initialVector.add(
				value.getNodeFunctional()
				     .evaluate(initialDisplacementValues()), nEulerian + key));
		lagrangian
			.getShapeFunctionMap()
			.forEach((key, value) -> previousVector.add(
				value.getNodeFunctional()
				     .evaluate(initialDisplacementValues()), nEulerian + key));
		lagrangian
			.getShapeFunctionMap()
			.forEach((key, value) -> previousVector.add(
				-value.getNodeFunctional()
				      .evaluate(initialDisplacementDerivatives()) * dt,
				nEulerian + key));
	}
	
	public void writeBoundaryValues(final SparseMatrix currentMatrix,
	                                final DenseVector rightHandSide,
	                                final double t)
	{
		eulerian.writeBoundaryValuesTo(new ComposedMixedFunction(boundaryVelocityValues(t)),
		                               f -> f.center()
		                                     .x() == 0,
		                               (f, sf) -> sf.hasVelocityFunction(),
		                               currentMatrix,
		                               rightHandSide);
		final int firstPressure = eulerian
			.getShapeFunctionMap()
			.values()
			.stream()
			.filter(MixedFunction::hasPressureFunction)
			.mapToInt(QkQkFunction::getGlobalIndex)
			.findFirst()
			.orElse(0);
		currentMatrix.deleteRow(firstPressure);
		currentMatrix.set(1, firstPressure, firstPressure);
		rightHandSide.set(0, firstPressure);
	}
	
	public void writeABf(final SparseMatrix constantMatrix)
	{
		final double alpha = rhoF / dt;
		final MixedCellIntegral<TPCell, ContinuousTPShapeFunction, ContinuousTPVectorFunction, QkQkFunction>
			mass = MixedCellIntegral.fromVelocityIntegral(
			new TPVectorCellIntegral<>(ScalarFunction.constantFunction(alpha),
			                           TPVectorCellIntegral.VALUE_VALUE));
		final MixedCellIntegral<TPCell, ContinuousTPShapeFunction, ContinuousTPVectorFunction, QkQkFunction>
			symGrad = MixedCellIntegral.fromVelocityIntegral(
			new TPVectorCellIntegral<>(ScalarFunction.constantFunction(nu), TPVectorCellIntegral.SYM_GRAD));
		final MixedTPCellIntegral<ContinuousTPShapeFunction, ContinuousTPVectorFunction, QkQkFunction> divValue
			= new MixedTPCellIntegral<>(ScalarFunction.constantFunction(-1), MixedTPCellIntegral.DIV_VALUE);
		final SparseMatrix ABf = new SparseMatrix(nEulerian, nEulerian);
		eulerianAlphaMass = new SparseMatrix(nEulerian, nEulerian);
		eulerian.writeCellIntegralsToMatrix(List.of(mass, symGrad, divValue), ABf);
		eulerian.writeCellIntegralsToMatrix(List.of(mass), eulerianAlphaMass);
		constantMatrix.addSmallMatrixInPlaceAt(ABf, 0, 0);
	}
	
	public void writeAs(final SparseMatrix constantMatrix)
	{
		final double beta = (rhoS - rhoF) / (dt * dt);
		final TPVectorCellIntegral<ContinuousTPVectorFunction> mass =
			new TPVectorCellIntegral<>(ScalarFunction.constantFunction(beta),
			                           TPVectorCellIntegral.VALUE_VALUE);
		final TPVectorCellIntegral<ContinuousTPVectorFunction> elast =
			new TPVectorCellIntegral<>(ScalarFunction.constantFunction(kappa),
			                           elastMethod);
		final SparseMatrix As = new SparseMatrix(nLagrangian, nLagrangian);
		lagrangian.writeCellIntegralsToMatrix(List.of(mass, elast), As);
		lagrangianBetaMass = new SparseMatrix(nLagrangian, nLagrangian);
		lagrangian.writeCellIntegralsToMatrix(List.of(mass), lagrangianBetaMass);
		constantMatrix.addSmallMatrixInPlaceAt(As, nEulerian, nEulerian);
	}
	
	public void writeCf(final SparseMatrix currentMatrix, final DenseVector currentIterate)
	{
		final SparseMatrix Cf = new SparseMatrix(nTransfer, nEulerian);
		int i = 0;
		for (final Map.Entry<Integer, QkQkFunction> sfEntry : eulerian.getShapeFunctionMap()
		                                                              .entrySet())
		{
			if (i++ % 40 == 0) System.out.println(100.0 * i / nEulerian);
			final QkQkFunction sf = sfEntry.getValue();
			final FunctionOnCells vOfX = concatenateVelocityWithX(sf, currentIterate);
			
			final TPVectorRightHandSideIntegral<ContinuousTPVectorFunction> transferEulerian
				= new TPVectorRightHandSideIntegral<>(vOfX,
				                                      TPVectorRightHandSideIntegral.H1);
			final DenseVector column = new DenseVector(nTransfer);
			if (sf.hasVelocityFunction())
			{
				lagrangian.writeCellIntegralsToRhs(List.of(transferEulerian), column, (K, lambd) ->
					                                   sf.getVelocityFunction()
					                                     .getNodeFunctionalPoint()
					                                     .sub(K.center())
					                                     .euclidianNorm() < 1. / nEulerCells + 1. / nLagrangianCells
				                                  );
				Cf.addColumn(column.mul(1), sfEntry.getKey());
			}
		}
		currentMatrix.addSmallMatrixInPlaceAt(Cf.mul(1), nEulerian + nLagrangian, 0);
		currentMatrix.addSmallMatrixInPlaceAt(Cf.transpose()
		                                        .mul(1), 0, nEulerian + nLagrangian);
	}
	
	private VectorTPFESpaceFunction<ContinuousTPVectorFunction> getX(final DenseVector currentIterate)
	{
		return new VectorTPFESpaceFunction<>(lagrangian.getShapeFunctionMap(),
		                                     getLagrangianIterate(currentIterate));
	}
	
	private MixedTPFESpaceFunction<QkQkFunction> getUp(final DenseVector currentIterate)
	{
		return new MixedTPFESpaceFunction<>(eulerian.getShapeFunctionMap(),
		                                    getEulerianfIterate(currentIterate));
	}
	
	private VectorFunctionOnCells<TPCell, TPFace> concatenateVelocityWithX(final MixedFunction sf,
	                                                                       final DenseVector currentIterate)
	{
		return VectorFunctionOnCells.concatenate(sf.getVelocityFunction(), getX(currentIterate));
//		final VectorTPFESpaceFunction<ContinuousTPVectorFunction> X = getX(currentIterate);
//		return new VectorFunctionOnCells<TPCell, TPFace>()
//		{
//			@Override
//			public int getRangeDimension()
//			{
//				return 2;
//			}
//
//			@Override
//			public CoordinateVector valueInCell(final CoordinateVector pos, final TPCell cell)
//			{
//				return sf.getVelocityFunction()
//				         .value(X.valueInCell(pos, cell));
//			}
//
//			@Override
//			public CoordinateMatrix gradientInCell(final CoordinateVector pos, final TPCell cell)
//			{
//
//				return sf.getVelocityFunction()
//				         .gradient(X.value(pos))
//				         .mmMul(X.gradientInCell(pos,
//				                                 cell));
//			}
//
//			@Override
//			public int getDomainDimension()
//			{
//				return 2;
//			}
//
//			@Override
//			public CoordinateVector value(final CoordinateVector pos)
//			{
//				return sf.getVelocityFunction()
//				         .value(X.value(pos));
//			}
//
//			@Override
//			public CoordinateMatrix gradient(final CoordinateVector pos)
//			{
//				return sf.getVelocityFunction()
//				         .gradient(X.value(pos))
//				         .mmMul(X.gradient(pos));
//			}
//		};
	}
	
	private DenseVector getLagrangianIterate(final DenseVector currentIterate)
	{
		if (currentIterate.getLength() != nEulerian + nLagrangian + nTransfer)
			throw new IllegalArgumentException("give currentIterate");
		return currentIterate.slice(nEulerian, nEulerian + nLagrangian);
	}
	
	private DenseVector getEulerianfIterate(final DenseVector currentIterate)
	{
		if (currentIterate.getLength() != nEulerian + nLagrangian + nTransfer)
			throw new IllegalArgumentException("give currentIterate");
		return currentIterate.slice(0, nEulerian);
	}
	
	private DenseVector getTransferIterate(final DenseVector currentIterate)
	{
		if (currentIterate.getLength() != nEulerian + nLagrangian + nTransfer)
			throw new IllegalArgumentException("give currentIterate");
		return currentIterate.slice(nEulerian + nLagrangian, nEulerian + nLagrangian + nTransfer);
	}
	
	public void writeCs(final SparseMatrix constantMatrix)
	{
		final SparseMatrix Cs = new SparseMatrix(nTransfer, nLagrangian);
		final TPVectorCellIntegral<ContinuousTPVectorFunction> transferLagrangian =
			new TPVectorCellIntegral<>(ScalarFunction.constantFunction(-1),
			                           TPVectorCellIntegral.H1);
		lagrangian.writeCellIntegralsToMatrix(List.of(transferLagrangian), Cs);
		lagrangianDual = Cs;
		constantMatrix.addSmallMatrixInPlaceAt(Cs.mul(1. / dt), nEulerian + nLagrangian, nEulerian);
		constantMatrix.addSmallMatrixInPlaceAt(Cs.transpose(), nEulerian, nEulerian + nLagrangian);
	}
	
	public void writeF(final DenseVector currentVector, final DenseVector currentIterate)
	{
		final DenseVector f = eulerianAlphaMass.mvMul(getEulerianfIterate(currentIterate));
		currentVector.addSmallVectorAt(f, 0);
	}
	
	public void writeG(final DenseVector currentVector,
	                   final DenseVector currentIterate,
	                   final DenseVector lastIterate)
	{
		final DenseVector g = lagrangianBetaMass.mvMul(
			getLagrangianIterate(currentIterate).mul(2)
			                                    .sub(getLagrangianIterate(lastIterate)));
		currentVector.addSmallVectorAt(g, nEulerian);
	}
	
	public void writeD(final DenseVector currentVector, final DenseVector currentIterate)
	{
		final DenseVector d = lagrangianDual.mvMul(getLagrangianIterate(currentIterate))
		                                    .mul(1. / dt);
		currentVector.addSmallVectorAt(d, nEulerian + nLagrangian);
	}
}
