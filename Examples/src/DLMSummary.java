import basic.PlotWindow;
import basic.ScalarFunction;
import basic.VectorFunction;
import distorted.*;
import distorted.geometry.DistortedCell;
import linalg.*;
import mixed.*;
import tensorproduct.*;
import tensorproduct.geometry.TPCell;

import java.util.List;
import java.util.Map;

public class DLMSummary
{
	final double dt;
	final int timeSteps;
	final double rhoF;
	final double rhoS;
	final double nu;
	final double kappa;
	TaylorHoodSpace eulerian;
	DistortedVectorSpace lagrangian;
	SparseMatrix lagrangianBetaMass;
	SparseMatrix eulerianAlphaMass;
	int nEulerian;
	int nLagrangian;
	int nTransfer;
	int eulerianPointsPerDimension = 40;
	List<CoordinateVector> eulerianPoints;
	
	public DLMSummary()
	{
		rhoF = 1;
		rhoS = 2;
		nu = 10;
		kappa = 10;
		dt = 0.03;
		timeSteps = 2;
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
		System.out.println(lastIterate.slice(nEulerian, nEulerian + nLagrangian));
		System.out.println(currentIterate.slice(nEulerian, nEulerian + nLagrangian));
		System.out.println("firstIts");
		writeBoundaryValues(constantSystemMatrix, currentIterate, 0);
		writeBoundaryValues(constantSystemMatrix, lastIterate, -dt);
		System.out.println("bdrs");
		final Map<CoordinateVector, Double> pvals = (new MixedFESpaceFunction<>(eulerian.getShapeFunctions(),
		                                                                        lastIterate).pressureValuesInPointsAtTime(
			eulerianPoints, -dt));
		final Map<CoordinateVector, CoordinateVector> vvals = (new MixedFESpaceFunction<>(
			eulerian.getShapeFunctions(), lastIterate).velocityValuesInPointsAtTime(eulerianPoints, -dt));
		
		final PlotWindow p = new PlotWindow();
		p.addPlot(new MatrixPlot(constantSystemMatrix));
		p.addPlot(new MatrixPlot(constantSystemMatrix.slice(new IntCoordinates(nEulerian, nEulerian),
		                                                    new IntCoordinates(nEulerian + nLagrangian,
		                                                                       nEulerian + nLagrangian))));
		for (int i = 1; i < timeSteps; i++)
		{
			rightHandSide = new DenseVector(constantRightHandSide);
			systemMatrix = new SparseMatrix(constantSystemMatrix);
			pvals.putAll(new MixedFESpaceFunction<>(eulerian.getShapeFunctions(),
			                                        currentIterate).pressureValuesInPointsAtTime(
				eulerianPoints, (i - 1) * dt));
			vvals.putAll(new MixedFESpaceFunction<>(eulerian.getShapeFunctions(),
			                                        currentIterate).velocityValuesInPointsAtTime(
				eulerianPoints, (i - 1) * dt));
			System.out.println("saved Iterate");
			writeF(rightHandSide, currentIterate.slice(0, nEulerian));
			System.out.println("f");
			System.out.println("RHS beforeG" + rightHandSide.slice(nEulerian, nEulerian + nLagrangian));
			writeG(rightHandSide, currentIterate.slice(nEulerian, nEulerian + nLagrangian),
			       lastIterate.slice(nEulerian, nEulerian + nLagrangian));
			System.out.println("RHS afterG" + rightHandSide.slice(nEulerian, nEulerian + nLagrangian));
			System.out.println("g");
			writeD(rightHandSide, currentIterate.slice(nEulerian, nEulerian + nLagrangian));
			System.out.println("d");
			writeCf(systemMatrix, currentIterate.slice(nEulerian, nEulerian + nLagrangian));
			System.out.println("cf");
			writeBoundaryValues(systemMatrix, currentIterate, i * dt);
			System.out.println("bdr");
			System.out.println(i + "th iteration");
			lastIterate = new DenseVector(currentIterate);
			currentIterate = systemMatrix.solve(rightHandSide);
			System.out.println("solved");
			System.out.println(
				"solved current iterate" + currentIterate.slice(nEulerian, nEulerian + nLagrangian));
			if (i ==1)
			{
				p.addPlot(new MatrixPlot(systemMatrix));
				p.addPlot(new MatrixPlot(systemMatrix.slice(new IntCoordinates(nEulerian, nEulerian),
				                                            new IntCoordinates(nEulerian + nLagrangian,
				                                                               nEulerian + nLagrangian))));
			}
		}
		
		final MixedPlot2DTime plot = new MixedPlot2DTime(pvals, vvals, eulerianPointsPerDimension);
		p.addPlot(plot.addOverlay(new CircleOverlay(lagrangian, plot)));
		p.addPlot(plot.addOverlay(new CircleOverlay(lagrangian, plot)));
	}
	
	public static void main(final String[] args)
	{
		final DLMSummary summary = new DLMSummary();
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
				return new CoordinateVector(2);
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
				return CoordinateVector.fromValues(pos.y(), 0);
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
				return new CoordinateVector(2);
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
				return CoordinateVector.fromValues(-7000, -7000);
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
		final IntCoordinates cellCounts = new IntCoordinates(3, 3);
		eulerian = new TaylorHoodSpace(startCoordinates, endCoordinates, cellCounts);
		eulerian.assembleCells();
		eulerian.assembleFunctions(1);
		nEulerian = eulerian.getShapeFunctions().size();
		eulerianPoints = eulerian.generatePlotPoints(eulerianPointsPerDimension);
	}
	
	public void initializeLagrangian()
	{
		final CoordinateVector center = CoordinateVector.fromValues(0.5, 0.5);
		final double radius = 0.2;
		lagrangian = new DistortedVectorSpace(center, radius, 1);
		lagrangian.assembleCells();
		lagrangian.assembleFunctions(1);
		nLagrangian = lagrangian.getShapeFunctions().size();
		nTransfer = nLagrangian;
	}
	
	public void generateFirstEulerianIterates(final DenseVector initialVector)
	{
		final MixedRightHandSideIntegral<TPCell, ContinuousTPShapeFunction, ContinuousTPVectorFunction, QkQkFunction>
			initialPressure = MixedRightHandSideIntegral.fromPressureIntegral(
			new TPRightHandSideIntegral<>(initialPressureValues(), TPRightHandSideIntegral.VALUE));
		final MixedRightHandSideIntegral<TPCell, ContinuousTPShapeFunction, ContinuousTPVectorFunction, QkQkFunction>
			initialVelocity = MixedRightHandSideIntegral.fromVelocityIntegral(
			new TPVectorRightHandSideIntegral<>(initialVelocityValues(),
			                                    TPVectorRightHandSideIntegral.VALUE));
		final DenseVector discreteInitialEulerian = new DenseVector(nEulerian);
		eulerian.writeCellIntegralsToRhs(List.of(initialPressure, initialVelocity), discreteInitialEulerian);
		initialVector.addSmallVectorAt(discreteInitialEulerian, 0);
	}
	
	public void generateFirstLagrangianIterates(final DenseVector initialVector, final DenseVector previousVector)
	{
		final DistortedVectorRightHandSideIntegral initialLagrangian = new DistortedVectorRightHandSideIntegral(
			initialDisplacementValues(), DistortedVectorRightHandSideIntegral.VALUE);
		final DistortedVectorRightHandSideIntegral initialLagrangianDerivative
			= new DistortedVectorRightHandSideIntegral(initialDisplacementDerivatives(),
			                                           DistortedVectorRightHandSideIntegral.VALUE);
		final DenseVector discreteInitialLagrangian = new DenseVector(nLagrangian);
		lagrangian.writeCellIntegralsToRhs(List.of(initialLagrangian), discreteInitialLagrangian);
		final DenseVector discreteLagrangianDerivative = new DenseVector(nLagrangian);
		lagrangian.writeCellIntegralsToRhs(List.of(initialLagrangianDerivative), discreteLagrangianDerivative);
		initialVector.addSmallVectorAt(discreteInitialLagrangian, nEulerian);
		final DenseVector previousLagrangian = discreteInitialLagrangian.sub(
			discreteLagrangianDerivative.mul(dt));
		previousVector.addSmallVectorAt(previousLagrangian, nEulerian);
	}
	
	public void writeBoundaryValues(final SparseMatrix currentMatrix, final DenseVector currentIterate, final double t)
	{
		eulerian.writeBoundaryValuesTo(new MixedFunction(boundaryVelocityValues(t)), currentMatrix,
		                               currentIterate);
		final int firstPressure = eulerian
			.getShapeFunctions()
			.values()
			.stream()
			.filter(MixedFunction::hasPressureFunction)
			.mapToInt(QkQkFunction::getGlobalIndex)
			.findFirst()
			.orElse(0);
		currentMatrix.deleteRow(firstPressure);
		currentMatrix.set(1, firstPressure, firstPressure);
		currentIterate.set(0, firstPressure);
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
		constantMatrix.addSmallMatrixAt(ABf, 0, 0);
	}
	
	public void writeAs(final SparseMatrix constantMatrix)
	{
		final double beta = (rhoS - rhoF) / dt;
		final double gamma = kappa * dt;
		final DistortedVectorCellIntegral mass = new DistortedVectorCellIntegral(beta,
		                                                                         DistortedVectorCellIntegral.VALUE_VALUE);
		final DistortedVectorCellIntegral gradGrad = new DistortedVectorCellIntegral(gamma,
		                                                                             DistortedVectorCellIntegral.GRAD_GRAD);
		final SparseMatrix As = new SparseMatrix(nLagrangian, nLagrangian);
		lagrangian.writeCellIntegralsToMatrix(List.of(mass, gradGrad), As);
		lagrangianBetaMass = new SparseMatrix(nLagrangian, nLagrangian);
		lagrangian.writeCellIntegralsToMatrix(List.of(mass), lagrangianBetaMass);
		constantMatrix.addSmallMatrixAt(As, nEulerian, nEulerian);
	}
	
	public void writeCf(final SparseMatrix currentMatrix, final Vector lagrangianIterate)
	{
		final DistortedVectorFESpaceFunction X = new DistortedVectorFESpaceFunction(
			lagrangian.getShapeFunctions(), lagrangianIterate);
		final SparseMatrix Cf = new SparseMatrix(nTransfer, nEulerian);
		int i = 0;
		for (final Map.Entry<Integer, QkQkFunction> sfEntry : eulerian.getShapeFunctions().entrySet())
		{
			if (i++ % 40 == 0) System.out.println(100.0 * i / nEulerian);
			final QkQkFunction sf = sfEntry.getValue();
			final DistortedVectorFunction vOfX = new DistortedVectorFunction()
			{
				
				@Override
				public CoordinateVector valueOnReferenceCell(final CoordinateVector pos, final DistortedCell cell)
				{
					return sf.getVelocityFunction().value(X.valueOnReferenceCell(pos, cell));
				}
				
				@Override
				public CoordinateMatrix gradientOnReferenceCell(final CoordinateVector pos, final DistortedCell cell)
				{
					return sf
						.getVelocityShapeFunction()
						.gradient(X.valueOnReferenceCell(pos, cell))
						.mmMul(X.gradientOnReferenceCell(pos, cell));
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
					throw new UnsupportedOperationException("Should not be evaluated");
				}
			};
			
			final DistortedFESpaceVectorRightHandSideIntegral transferEulerian
				= new DistortedFESpaceVectorRightHandSideIntegral(vOfX,
				                                                  DistortedRightHandSideIntegral.H1);
			final DenseVector column = new DenseVector(nTransfer);
			if (sf.hasVelocityFunction())
			
			{
				lagrangian.writeCellIntegralsToRhs(List.of(transferEulerian), column);
				Cf.addColumn(column, sfEntry.getKey());
			}
		}
		currentMatrix.addSmallMatrixAt(Cf, nEulerian + nLagrangian, 0);
		currentMatrix.addSmallMatrixAt(Cf.transpose(), 0, nEulerian + nLagrangian);
	}
	
	public void writeCs(final SparseMatrix constantMatrix)
	{
		final SparseMatrix Cs = new SparseMatrix(nTransfer, nLagrangian);
		final DistortedVectorCellIntegral transferLagrangian = new DistortedVectorCellIntegral(-1,
		                                                                                       DistortedVectorCellIntegral.H1);
		lagrangian.writeCellIntegralsToMatrix(List.of(transferLagrangian), Cs);
		constantMatrix.addSmallMatrixAt(Cs, nEulerian + nLagrangian, nEulerian);
		constantMatrix.addSmallMatrixAt(Cs.transpose(), nEulerian, nEulerian + nLagrangian);
	}
	
	public void writeF(final DenseVector currentVector, final Vector eulerianIterate)
	{
		final DenseVector f = eulerianAlphaMass.mvMul(eulerianIterate);
		currentVector.addSmallVectorAt(f, 0);
	}
	
	public void writeG(final DenseVector currentVector, final Vector lagrangianIterate, final Vector lastLagrangianIterate)
	{
		final DenseVector g = lagrangianBetaMass.mvMul(
			lagrangianIterate.mul(2).sub(lastLagrangianIterate).mul(1. / dt));
		currentVector.addSmallVectorAt(g, nEulerian);
	}
	
	public void writeD(final DenseVector currentVector, final Vector lagrangianIterate)
	{
		final DenseVector d = lagrangianBetaMass.mvMul(lagrangianIterate).mul(-1. / (rhoS - rhoF));
		currentVector.addSmallVectorAt(d, nEulerian + nLagrangian);
	}
}
