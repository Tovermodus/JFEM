import basic.*;
import distorted.*;
import distorted.geometry.DistortedCell;
import linalg.*;
import mixed.*;
import tensorproduct.ContinuousTPShapeFunction;
import tensorproduct.ContinuousTPVectorFunction;
import tensorproduct.TPVectorCellIntegral;
import tensorproduct.geometry.TPCell;

import java.util.List;
import java.util.Map;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.stream.Collectors;

public class DLMElast
{
	final double dt;
	final int timeSteps;
	final double rhoF;
	final double rhoS;
	final double nu;
	final double lameLamb;
	final double lameMu;
	TaylorHoodSpace eulerian;
	DistortedVectorSpace lagrangian;
	SparseMatrix lagrangianMass;
	SparseMatrix eulerianMass;
	SparseMatrix Cf;
	SparseMatrix lagrangianH1;
	int nEulerian;
	int nLagrangian;
	int nTransfer;
	int eulerianPointsPerDimension = 40;
	int nEulerCells = 14;
	int nLagrangeRefines = 1;
	List<CoordinateVector> eulerianPoints;
	private final int lagrangeDegree = 2;
	private final int eulerDegree = 1;
	private final String elastMethod = DistortedVectorCellIntegral.GRAD_GRAD;
	
	public DLMElast()
	{
		final PerformanceArguments.PerformanceArgumentBuilder builder =
			new PerformanceArguments.PerformanceArgumentBuilder();
		builder.executeChecks = false;
		builder.build();
		rhoF = 1;
		rhoS = 2;
		lameLamb = 10;
		lameMu = 1;
		nu = 1;
		dt = 0.01;
		timeSteps = 5;
		initializeEulerian();
		initializeLagrangian();
		
		final ExecutorService ex = Executors.newSingleThreadExecutor();
		final Interruptor interruptor = new Interruptor();
		ex.execute(interruptor);
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
		final DistortedVectorFESpaceFunction X = getX(lastIterate);
		final Map<CoordinateVector, Double> dvals = eulerianPoints
			.stream()
			.map(p ->
			     {
				     double val = 1;
				     CoordinateVector Xp = X.value(p);
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
		final Map<CoordinateVector, CoordinateVector> LfXVals = concatenateVelocityWithX(getUp(lastIterate),
		                                                                                 lastIterate)
			.valuesInPointsAtTime(eulerianPoints, -dt);
		
		final DenseMatrix iterateHistory = new DenseMatrix(timeSteps + 1, nEulerian + nLagrangian + nTransfer);
		final DenseMatrix rhsHistory = new DenseMatrix(timeSteps + 1, nEulerian + nLagrangian + nTransfer);
		iterateHistory.addRow(lastIterate, 0);
		iterateHistory.addRow(currentIterate, 1);
		final SparseMatrix precond = new SparseMatrix(constantSystemMatrix);
		//new DenseMatrix(precond).get
		writeCf(precond, currentIterate);
		final BlockSparseMatrix precondBlocks = new BlockSparseMatrix(precond, new int[]{0, nEulerian});
		final DenseMatrix lagInverse = new BlockDenseMatrix(new DenseMatrix(precondBlocks.getBlockMatrix(1,
		                                                                                                 1)),
		                                                    1).getInvertedDiagonalMatrix()
		                                                      .toDense();
		System.out.println("laginv");
		final DenseMatrix firstmmul = precondBlocks.getBlockMatrix(0, 1)
		                                           .mmMul(lagInverse);
		System.out.println("laginv");
		final DenseMatrix secondmmul = firstmmul
			.mmMul(precondBlocks.getBlockMatrix(1, 0));
		System.out.println("laginv");
		final BlockDenseMatrix schur = new BlockDenseMatrix(precondBlocks.getBlockMatrix(0, 0)
		                                                                 .add(secondmmul),
		                                                    4);
		final BlockDenseMatrix schurBlockInverse = schur.getInvertedDiagonalMatrix();
		//final DenseMatrix schurDense = schur.toDense();
		//final DenseMatrix schurDenseInverse = schurDense.inverse();
		//System.out.println("############");
		//final DenseVector r = new DenseVector(schurDense.getVectorSize());
//		System.out.println(schur.sub(schurDense)
//		                        .absMaxElement());
//		System.out.println(schurDenseInverse.mmMul(schurDense)
//		                                    .sub(DenseMatrix.identity(schur.getCols()))
//		                                    .absMaxElement());
//		System.out.println(schurDenseInverse.mmMul(schur)
//		                                    .sub(DenseMatrix.identity(schur.getCols()))
//		                                    .absMaxElement());
//		System.out.println("##########");
//		System.out.println(schurBlockInverse.getShape() + " " + schur.getShape());
//		System.out.println("mmul");
		final VectorMultiplyable precondInverse = new VectorMultiplyable()
		{
			final IterativeSolver it = new IterativeSolver();
			final SparseMvMul l = new SparseMvMul(precondBlocks.getBlockMatrix(0, 1));
			final SparseMvMul inv = new SparseMvMul(lagInverse);
			final SparseMvMul sch = new SparseMvMul(schur);
			final SparseMvMul schInv = new SparseMvMul(schurBlockInverse);
			final SparseMvMul r = new SparseMvMul(precondBlocks.getBlockMatrix(1, 0));
			
			@Override
			public int getVectorSize()
			{
				return nEulerian + nLagrangian + nTransfer;
			}
			
			@Override
			public int getTVectorSize()
			{
				throw new UnsupportedOperationException("not implemented yet");
			}
			
			@Override
			public Vector mvMul(final Vector vector)
			{
				final Vector g =
					vector.slice(0, nEulerian)
					      .sub(l.mvMul(inv.mvMul(vector.slice(nEulerian,
					                                          nEulerian + nLagrangian + nTransfer))));
				
				//it = new IterativeSolver();
				final DenseVector eul = //new DenseVector(schurDense.mvMul(g));
//					(DenseVector) it.solvePGMRES(sch, schInv, g,
//					                             1e-10);
					schInv.mvMul(g);
				System.out.println(sch.mvMul(eul)
				                      .sub(g)
				                      .absMaxElement() + " subCG1");
				final DenseVector lag = inv.mvMul(vector.slice(nEulerian,
				                                               nEulerian + nLagrangian + nTransfer)
				                                        .sub(r.mvMul(eul)));
				final DenseVector ret = new DenseVector(vector.getLength());
				ret.addSmallVectorAt(eul, 0);
				ret.addSmallVectorAt(lag, nEulerian);
				return ret;
			}
			
			@Override
			public Vector tvMul(final Vector vector)
			{
				throw new UnsupportedOperationException("not implemented yet");
			}
		};
		System.out.println("created preconditioner");
		int i;
		final IterativeSolver itp = new IterativeSolver(true);
		for (i = 1; i < timeSteps && interruptor.isRunning(); i++)
		{
			rightHandSide = new DenseVector(constantRightHandSide);
			systemMatrix = new SparseMatrix(constantSystemMatrix);
			final int drawInterval;
			if (timeSteps > 300)
				drawInterval = timeSteps / 300;
			else
				drawInterval = 1;
			if (i % drawInterval == 0 || i < 10)
				writeOutVals(currentIterate, lastIterate, pvals, vvals, dvals, derivVals, uXvals,
				             LfXVals, i);
			System.out.println("saved Iterate");
			writeF(rightHandSide, currentIterate);
			System.out.println("f");
			writeG(rightHandSide, currentIterate, lastIterate);
			System.out.println("g");
			writeD(rightHandSide, currentIterate);
			System.out.println("d");
			writeCf(systemMatrix, currentIterate);
			System.out.println("cf");
			writeAfSemiImplicit(systemMatrix, currentIterate);
			System.out.println("AfImplicit");
			writeBoundaryValues(systemMatrix, rightHandSide, i * dt);
			System.out.println("bdr");
			System.out.println(i + "th iteration");
			lastIterate = new DenseVector(currentIterate);
			rhsHistory.addRow(rightHandSide, i + 1);
			//final IterativeSolver it = new IterativeSolver();
			//it.showProgress = true;
			//it.solveGMRES(systemMatrix, rightHandSide, 1e-7);
			itp.showProgress = true;
			System.out.println(systemMatrix.sub(systemMatrix.transpose())
			                               .absMaxElement() + "ABBBBBB");
//			currentIterate = (DenseVector) itp.solvePGMRES(new SparseMvMul(systemMatrix),
//			                                               precondInverse,
//			                                               rightHandSide,
//			                                               1e-7);
			currentIterate = systemMatrix.solve(rightHandSide);
			//System.out.println(iterrr.sub(currentIterate).absMaxElement());
			//	System.out.println("newit" + getLagrangianIterate(currentIterate));
			iterateHistory.addRow(currentIterate, i + 1);
			System.out.println("solved");
		}
		interruptor.running = false;
		final PlotWindow p = new PlotWindow();
		p.addPlot(new MatrixPlot(iterateHistory, "iterateHistory"));
		p.addPlot(new MatrixPlot(rhsHistory, "rhsHistory"));
		final MixedPlot2DTime plot = new MixedPlot2DTime(pvals, vvals, eulerianPointsPerDimension, "VVals");
		final MixedPlot2DTime plot2 = new MixedPlot2DTime(pvals, uXvals, eulerianPointsPerDimension, "u at X");
		final MixedPlot2DTime plot4 = new MixedPlot2DTime(pvals, LfXVals, eulerianPointsPerDimension, "LfX");
		final MixedPlot2DTime plot3 = new MixedPlot2DTime(pvals, derivVals, eulerianPointsPerDimension,
		                                                  "derivative of X");
		//System.out.println("plotttt");
//		p.addPlot(new MixedPlot2DTime(pvals, uXvals, eulerianPointsPerDimension, "uXVals"));
		try
		{
			Thread.sleep(1000);
		} catch (final InterruptedException e)
		{
			e.printStackTrace();
		}
//		p.addPlot(new MixedPlot2DTime(pvals, derivVals, eulerianPointsPerDimension, "derivVals"));
		
		final ScalarPlot2DTime disPLacementplot = new ScalarPlot2DTime(dvals, eulerianPointsPerDimension,
		                                                               "dVals");
		p.addPlot(disPLacementplot);
		p.addPlot(plot);
		p.addPlot(plot2);
		p.addPlot(plot4);
		p.addPlot(plot3);
		plot.addOverlay(new DistortedOverlay(eulerianPoints, lagrangian,
		                                     iterateHistory.slice(new IntCoordinates(0, nEulerian),
		                                                          new IntCoordinates(i,
		                                                                             nEulerian + nLagrangian)),
		                                     7));
	}
	
	private void writeOutVals(final DenseVector currentIterate,
	                          final DenseVector lastIterate,
	                          final Map<CoordinateVector, Double> pvals,
	                          final Map<CoordinateVector, CoordinateVector> vvals,
	                          final Map<CoordinateVector, Double> dvals,
	                          final Map<CoordinateVector, CoordinateVector> derivVals,
	                          final Map<CoordinateVector, CoordinateVector> uXvals,
	                          final Map<CoordinateVector, CoordinateVector> LfXVals,
	                          final int i)
	{
		final DistortedVectorFESpaceFunction X;
		pvals.putAll(getUp(currentIterate).pressureValuesInPointsAtTime(eulerianPoints, (i - 1) * dt));
		vvals.putAll(getUp(currentIterate).velocityValuesInPointsAtTime(eulerianPoints, (i - 1) * dt));
		X = getX(currentIterate);
		final DistortedVectorFESpaceFunction finalX1 = X;
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
				                  return Map.entry(Xp.addCoordinate((i - 1) * dt), val);
			                  })
			             .collect(Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue)));
		final VectorFunction X_prime = new DistortedVectorFESpaceFunction(lagrangian.getShapeFunctions(),
		                                                                  getLagrangianIterate(
			                                                                  currentIterate.sub(
				                                                                  lastIterate))
			                                                                  .mul(1. / dt));
		derivVals.putAll(getX(currentIterate.sub(lastIterate)
		                                    .mul(1. / dt))
			                 .valuesInPointsAtTime(eulerianPoints, (i - 1) * dt));
		
		final VectorFunction U_of_X = new DistortedVectorFESpaceFunction(lagrangian.getShapeFunctions(),
		                                                                 Cf.mvMul(getEulerianfIterate(
			                                                                 currentIterate)));
		
		LfXVals.putAll(U_of_X.valuesInPointsAtTime(eulerianPoints, (i - 1) * dt));
		
		uXvals.putAll(concatenateVelocityWithX(getUp(currentIterate), currentIterate)
			              .valuesInPointsAtTime(eulerianPoints, (i - 1) * dt));
	}
	
	public static void main(final String[] args)
	{
		new DLMElast();
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
		nEulerian = eulerian.getShapeFunctions()
		                    .size();
		eulerianPoints = eulerian.generatePlotPoints(eulerianPointsPerDimension);
	}
	
	public void initializeLagrangian()
	{
		final CoordinateVector center = CoordinateVector.fromValues(0.5, 0.5);
		final double radius = 0.3;
		lagrangian = new DistortedVectorSpace(center, radius, nLagrangeRefines);
		lagrangian.assembleCells();
		lagrangian.assembleFunctions(lagrangeDegree);
		nLagrangian = lagrangian.getShapeFunctions()
		                        .size();
		nTransfer = nLagrangian;
	}
	
	public void generateFirstEulerianIterates(final DenseVector initialVector)
	{
		final MixedFunction initialVelocityFunction = new MixedFunction(initialVelocityValues());
		final MixedFunction initialPressureFunction = new MixedFunction(initialPressureValues());
		eulerian
			.getShapeFunctions()
			.forEach((key, value) -> initialVector.add(
				value.getNodeFunctional()
				     .evaluate(initialVelocityFunction), key));
		eulerian
			.getShapeFunctions()
			.forEach((key, value) -> initialVector.add(
				value.getNodeFunctional()
				     .evaluate(initialPressureFunction), key));
	}
	
	public void generateFirstLagrangianIterates(final DenseVector initialVector, final DenseVector previousVector)
	{
		lagrangian
			.getShapeFunctions()
			.forEach((key, value) -> initialVector.add(
				value.getNodeFunctional()
				     .evaluate(initialDisplacementValues()), nEulerian + key));
		lagrangian
			.getShapeFunctions()
			.forEach((key, value) -> previousVector.add(
				value.getNodeFunctional()
				     .evaluate(initialDisplacementValues()), nEulerian + key));
		lagrangian
			.getShapeFunctions()
			.forEach((key, value) -> previousVector.add(
				-value.getNodeFunctional()
				      .evaluate(initialDisplacementDerivatives()) * dt,
				nEulerian + key));
	}
	
	public void writeBoundaryValues(final SparseMatrix currentMatrix,
	                                final DenseVector rightHandSide,
	                                final double t)
	{
		eulerian.writeBoundaryValuesTo(new MixedFunction(boundaryVelocityValues(t)),
		                               f -> f.center()
		                                     .x() == 0,
		                               (f, sf) -> sf.hasVelocityFunction(),
		                               currentMatrix,
		                               rightHandSide);
		final int firstPressure = eulerian
			.getShapeFunctions()
			.values()
			.stream()
			.filter(MixedFunction::hasPressureFunction)
			.mapToInt(QkQkFunction::getGlobalIndex)
			.findFirst()
			.orElse(0);
		currentMatrix.deleteColumn(firstPressure);
		currentMatrix.deleteRow(firstPressure);
		currentMatrix.set(1, firstPressure, firstPressure);
		rightHandSide.set(0, firstPressure);
	}
	
	public void writeABf(final SparseMatrix constantMatrix)
	{
		final MixedCellIntegral<TPCell, ContinuousTPShapeFunction, ContinuousTPVectorFunction, QkQkFunction>
			mass = MixedCellIntegral.fromVelocityIntegral(
			new TPVectorCellIntegral<>(ScalarFunction.constantFunction(rhoF),
			                           TPVectorCellIntegral.VALUE_VALUE));
		final MixedCellIntegral<TPCell, ContinuousTPShapeFunction, ContinuousTPVectorFunction, QkQkFunction>
			symGrad = MixedCellIntegral.fromVelocityIntegral(
			new TPVectorCellIntegral<>(ScalarFunction.constantFunction(nu), TPVectorCellIntegral.SYM_GRAD));
		final MixedTPCellIntegral<ContinuousTPShapeFunction, ContinuousTPVectorFunction, QkQkFunction> divValue
			= new MixedTPCellIntegral<>(ScalarFunction.constantFunction(-1), MixedTPCellIntegral.DIV_VALUE);
		final SparseMatrix ABf = new SparseMatrix(nEulerian, nEulerian);
		eulerianMass = new SparseMatrix(nEulerian, nEulerian);
		eulerian.writeCellIntegralsToMatrix(List.of(mass), ABf);
		eulerian.writeCellIntegralsToMatrix(List.of(mass), eulerianMass);
		ABf.mulInPlace(1. / dt);
		eulerian.writeCellIntegralsToMatrix(List.of(symGrad, divValue), ABf);
		constantMatrix.addSmallMatrixAt(ABf, 0, 0);
	}
	
	public void writeAfSemiImplicit(final SparseMatrix currentMatrix, final DenseVector currentIterate)
	{
//		final VectorFunctionOnCells<TPCell, TPFace> vel = getUp(currentIterate).getVelocityFunction();
//		final VectorFunction semiImplicitWeight1 = VectorFunction.fromLambda((x) -> vel.value(x)
//		                                                                               .mul(rhoF / 2), 2, 2);
//		final VectorFunction semiImplicitWeight2 = VectorFunction.fromLambda((x) -> vel.value(x)
//		                                                                               .mul(-rhoF / 2), 2, 2);
//		final TPVectorCellIntegralOnCell<ContinuousTPVectorFunction> convection1 =
//			new TPVectorCellIntegralOnCell<>(semiImplicitWeight1, TPVectorCellIntegral.GRAD_VALUE);
//		final TPVectorCellIntegralOnCell<ContinuousTPVectorFunction> convection2 =
//			new TPVectorCellIntegralOnCell<>(semiImplicitWeight2, TPVectorCellIntegral.VALUE_GRAD);
//		final MixedCellIntegral<TPCell, ContinuousTPShapeFunction, ContinuousTPVectorFunction,
//			QkQkFunction> mixedConvection1 = MixedCellIntegral.fromVelocityIntegral(convection1);
//		final MixedCellIntegral<TPCell, ContinuousTPShapeFunction, ContinuousTPVectorFunction,
//			QkQkFunction> mixedConvection2 = MixedCellIntegral.fromVelocityIntegral(convection2);
//		final SparseMatrix convection = new SparseMatrix(nEulerian, nEulerian);
//		eulerian.writeCellIntegralsToMatrix(List.of(mixedConvection1, mixedConvection2), convection);
//		currentMatrix.addSmallMatrixAt(convection, 0, 0);
	}
	
	public void writeAs(final SparseMatrix constantMatrix)
	{
		final DistortedVectorCellIntegral mass = new DistortedVectorCellIntegral(rhoS,
		                                                                         DistortedVectorCellIntegral.VALUE_VALUE);
		final DistortedVectorCellIntegral elastLamb = new DistortedVectorCellIntegral(lameLamb,
		                                                                              DistortedVectorCellIntegral.TRACE_SYM);
		final DistortedVectorCellIntegral elastMu = new DistortedVectorCellIntegral(lameMu,
		                                                                            DistortedVectorCellIntegral.SYM_SYM);
		final SparseMatrix As = new SparseMatrix(nLagrangian, nLagrangian);
		lagrangian.writeCellIntegralsToMatrix(List.of(mass), As);
		As.mulInPlace(1. / (dt * dt));
		lagrangian.writeCellIntegralsToMatrix(List.of(elastLamb, elastMu), As);
		lagrangianMass = new SparseMatrix(nLagrangian, nLagrangian);
		lagrangian.writeCellIntegralsToMatrix(List.of(mass), lagrangianMass);
//		System.out.println(As.transpose()
//		                     .sub(As)
//		                     .absMaxElement() + "AsAbsMax");
		constantMatrix.addSmallMatrixAt(As, nEulerian, nEulerian);
	}
	
	public void writeCf(final SparseMatrix currentMatrix, final DenseVector currentIterate)
	{
		Cf = new SparseMatrix(nTransfer, nEulerian);
		eulerian.getShapeFunctions()
		        .entrySet()
		        .stream()
		        .parallel()
		        .filter(e -> e.getValue()
		                      .hasVelocityFunction())
		        .forEach(
			        sfEntry ->
			        {
				        final QkQkFunction sf = sfEntry.getValue();
				        final DistortedVectorFunction vOfX = concatenateVelocityWithX(sf,
				                                                                      currentIterate);
				
				        final DistortedVectorDistortedRightHandSideIntegral transferEulerian
					        = new DistortedVectorDistortedRightHandSideIntegral(vOfX,
					                                                            DistortedVectorDistortedRightHandSideIntegral.VALUE);
				        final DenseVector column = new DenseVector(nTransfer);
				        lagrangian.writeCellIntegralsToRhs(List.of(transferEulerian), column,
				                                           (K, lambd) ->
					                                           sf.getVelocityShapeFunction()
					                                             .getNodeFunctionalPoint()
					                                             .sub(K.center())
					                                             .euclidianNorm() < 1 + 2. / nEulerCells + lagrangian
						                                           .getMaxDiam());
				        Cf.addColumn(column, sfEntry.getKey());
			        });
		currentMatrix.addSmallMatrixAt(Cf, nEulerian + nLagrangian, 0);
		currentMatrix.addSmallMatrixAt(Cf.transpose(), 0, nEulerian + nLagrangian);
	}
	
	private DistortedVectorFESpaceFunction getX(final DenseVector currentIterate)
	{
		return new DistortedVectorFESpaceFunction(lagrangian.getShapeFunctions(),
		                                          getLagrangianIterate(currentIterate));
	}
	
	private MixedTPFESpaceFunction<QkQkFunction> getUp(final DenseVector currentIterate)
	{
		return new MixedTPFESpaceFunction<>(eulerian.getShapeFunctions(), getEulerianfIterate(currentIterate));
	}
	
	private DistortedVectorFunction concatenateVelocityWithX(final MixedFunction sf,
	                                                         final DenseVector currentIterate)
	{
		
		final DistortedVectorFESpaceFunction X = getX(currentIterate);
		
		return new DistortedVectorFunction()
		{
			@Override
			public CoordinateVector valueOnReferenceCell(final CoordinateVector pos,
			                                             final DistortedCell cell)
			{
				return sf.getVelocityFunction()
				         .value(X.valueOnReferenceCell(pos, cell));
			}
			
			@Override
			public CoordinateMatrix gradientOnReferenceCell(final CoordinateVector pos,
			                                                final DistortedCell cell)
			{
				throw new IllegalArgumentException("askljhgd");
//				return sf
//					.getVelocityFunction()
//					.gradient(X.valueOnReferenceCell(pos, cell))
//					.mmMul(X.gradientOnReferenceCell(pos, cell));
			}
			
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
				return sf.getVelocityFunction()
				         .value(X.value(pos));
			}
		};
	}
	
	private DenseVector getLagrangianIterate(final Vector currentIterate)
	{
		if (currentIterate.getLength() != nEulerian + nLagrangian + nTransfer)
			throw new IllegalArgumentException("give currentIterate");
		return currentIterate.slice(nEulerian, nEulerian + nLagrangian);
	}
	
	private DenseVector getEulerianfIterate(final Vector currentIterate)
	{
		if (currentIterate.getLength() != nEulerian + nLagrangian + nTransfer)
			throw new IllegalArgumentException("give currentIterate");
		return currentIterate.slice(0, nEulerian);
	}
	
	private DenseVector getTransferIterate(final Vector currentIterate)
	{
		if (currentIterate.getLength() != nEulerian + nLagrangian + nTransfer)
			throw new IllegalArgumentException("give currentIterate");
		return currentIterate.slice(nEulerian + nLagrangian, nEulerian + nLagrangian + nTransfer);
	}
	
	public void writeCs(final SparseMatrix constantMatrix)
	{
		final SparseMatrix Cs = new SparseMatrix(nTransfer, nLagrangian);
		final DistortedVectorCellIntegral transferLagrangian
			= new DistortedVectorCellIntegral(DistortedVectorCellIntegral.VALUE_VALUE);
		lagrangian.writeCellIntegralsToMatrix(List.of(transferLagrangian), Cs);
		lagrangianH1 = Cs;
		constantMatrix.addSmallMatrixAt(Cs.mul(-1. / dt), nEulerian + nLagrangian, nEulerian);
		constantMatrix.addSmallMatrixAt(Cs.transpose()
		                                  .mul(-1), nEulerian, nEulerian + nLagrangian);
	}
	
	public void writeF(final DenseVector currentVector, final DenseVector currentIterate)
	{
		final DenseVector f = eulerianMass.mvMul(getEulerianfIterate(currentIterate))
		                                  .mul(1. / dt);
		currentVector.addSmallVectorAt(f, 0);
	}
	
	public void writeG(final DenseVector currentVector,
	                   final DenseVector currentIterate,
	                   final DenseVector lastIterate)
	{
		final DenseVector g = lagrangianMass.mvMul(
			getLagrangianIterate(currentIterate).mul(2)
			                                    .sub(getLagrangianIterate(lastIterate)));
		currentVector.addSmallVectorAt(g.mul(1. / (dt * dt)), nEulerian);
	}
	
	public void writeD(final DenseVector currentVector, final DenseVector currentIterate)
	{
		final DenseVector d = lagrangianH1.mvMul(getLagrangianIterate(currentIterate))
		                                  .mul(-1. / dt);
		currentVector.addSmallVectorAt(d, nEulerian + nLagrangian);
	}
	
	private void subtractIdentityOnLagrangianMesh(final CoordinateVector k, final CoordinateVector v)
	{
		final CoordinateVector keyWithoutTime = CoordinateVector.fromValues(k.x(), k.y());
		if (keyWithoutTime.sub(lagrangian.center)
		                  .euclidianNorm() < lagrangian.radius - 1e-2)
			v.subInPlace(keyWithoutTime);
	}
}
