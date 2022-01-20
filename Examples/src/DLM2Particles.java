import basic.CellIntegral;
import basic.PlotWindow;
import basic.RightHandSideIntegral;
import basic.ScalarFunction;
import distorted.*;
import distorted.geometry.DistortedCell;
import dlm.HyperbolicCartesianDistorted;
import linalg.*;
import mixed.*;
import scala.Function2;
import tensorproduct.CartesianGridSpace;
import tensorproduct.ContinuousTPShapeFunction;
import tensorproduct.ContinuousTPVectorFunction;
import tensorproduct.TPVectorCellIntegral;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentSkipListMap;
import java.util.function.Function;
import java.util.function.ToDoubleFunction;

public class DLM2Particles
	extends HyperbolicCartesianDistorted<QkQkFunction, ContinuousTPShapeFunction, ContinuousTPVectorFunction>
{
	
	final Map<CoordinateVector, CoordinateVector> velocityValues;
	final Map<CoordinateVector, Double> pressureValues;
	final Map<CoordinateVector, CoordinateVector> XValues;
	final Map<CoordinateVector, CoordinateVector> XPrimeValues;
	final Map<CoordinateVector, CoordinateVector> UXValues;
	final List<CoordinateVector> plotPoints;
	VectorMultiplyable precond;
	
	public DLM2Particles(final double dt,
	                     final int timeSteps,
	                     final CartesianGridSpace<QkQkFunction, MixedValue, MixedGradient, MixedHessian> backgroundSpace,
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
//		final IterativeSolver it = new IterativeSolver();
//		it.showProgress = true;
//		it.showInterrupt = true;
		return (A, b) ->
		{
			System.out.println(A.getClass());
			System.out.println(A.getRows());
			//final DenseMatrix Ad = new DenseMatrix(A);
			System.out.println(A.getRows());
			return ((SparseMatrix) A).solve(b);//new DenseVector(it.solvePGMRES(A, precond, b, 1e-6));
		};
	}
	
	VectorMultiplyable getPrecond()
	{
		final Matrix schur = getSchurComplement();
		final BlockDenseMatrix schurBlockInverse = new BlockDenseMatrix(schur, 4).getInvertedDiagonalMatrix();
		final Matrix B = getAllBackgroundTransferMatrices(0);
		final Matrix C = getAllBackgroundTransferTransposeMatrices(0);
		final Matrix lagInverse = new DenseMatrix(getAllParticleMatrices(0)).inverse();
		return new VectorMultiplyable()
		{
			final IterativeSolver it = new IterativeSolver();
			final SparseMvMul b = new SparseMvMul(B);
			final SparseMvMul inv = new SparseMvMul(lagInverse);
			final SparseMvMul sch = new SparseMvMul(schur);
			final SparseMvMul schInv = new SparseMvMul(schurBlockInverse);
			final SparseMvMul c = new SparseMvMul(C);
			
			@Override
			public int getVectorSize()
			{
				return getSystemSize();
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
					getBackGroundIterate(vector)
						.sub(b.mvMul(inv.mvMul(getAllParticleIterate(vector, 0))));
				
				//it = new IterativeSolver();
				final Vector eul = it.solveGMRES(sch, g, 1e-5);
				System.out.println(sch.mvMul(eul)
				                      .sub(g)
				                      .absMaxElement() + " subCG1");
				final DenseVector lag = inv.mvMul(getAllParticleIterate(vector, 0)
					                                  .sub(c.mvMul(eul)));
				final DenseVector ret = new DenseVector(vector.getLength());
				ret.addSmallVectorAt(eul, 0);
				ret.addSmallVectorAt(lag, (int) eul.size());
				return ret;
			}
			
			@Override
			public Vector tvMul(final Vector vector)
			{
				throw new UnsupportedOperationException("not implemented yet");
			}
		};
	}
	
	@Override
	protected void postIterationCallback()
	{
//		if (precond == null)
//		{
//			precond = getPrecond();
//		}
		final MixedFESpaceFunction<QkQkFunction, TPCell, TPFace> Up =
			new MixedTPFESpaceFunction<>(backgroundSpace.getShapeFunctions(),
			                             getBackGroundIterate(getCurrentIterate()));
		velocityValues.putAll(Up.velocityValuesInPointsAtTime(plotPoints, getTime()));
		pressureValues.putAll(Up.pressureValuesInPointsAtTime(plotPoints, getTime()));
		final DistortedVectorFESpaceFunction X =
			new DistortedVectorFESpaceFunction(particleSpaces.get(0)
			                                                 .getShapeFunctions(),
			                                   getParticleIterate(getCurrentIterate(), 0));
		XValues.putAll(X.valuesInPointsAtTime(plotPoints, getTime()));
		final DistortedVectorFESpaceFunction XPrime =
			new DistortedVectorFESpaceFunction(particleSpaces.get(0)
			                                                 .getShapeFunctions(),
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
			new TPVectorCellIntegral<ContinuousTPVectorFunction>(1, TPVectorCellIntegral.SYM_GRAD);
		final MixedCellIntegral<TPCell, ContinuousTPShapeFunction, ContinuousTPVectorFunction, QkQkFunction>
			symGradMixed = MixedTPCellIntegral.fromVelocityIntegral(symGrad);
		final MixedCellIntegral<TPCell, ContinuousTPShapeFunction, ContinuousTPVectorFunction, QkQkFunction>
			divValue
			= new MixedTPCellIntegral<>(ScalarFunction.constantFunction(-1), MixedTPCellIntegral.DIV_VALUE);
		return List.of(symGradMixed, divValue);
	}
	
	@Override
	protected List<CellIntegral<DistortedCell, DistortedVectorShapeFunction>> getParticleIntegrals(final int particleId)
	{
		final double lameLambda = 1;
		final double lameMu = 10;
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
			new TPVectorCellIntegral<ContinuousTPVectorFunction>(1, TPVectorCellIntegral.VALUE_VALUE);
		final MixedCellIntegral<TPCell, ContinuousTPShapeFunction, ContinuousTPVectorFunction, QkQkFunction>
			mixedMass
			= MixedCellIntegral.fromVelocityIntegral(mass);
		return List.of(mixedMass);
	}
	
	@Override
	protected List<CellIntegral<DistortedCell, DistortedVectorShapeFunction>> getParticleMassIntegrals(
		final int particleId)
	{
		final DistortedVectorCellIntegral mass = new DistortedVectorCellIntegral(100,
		                                                                         DistortedVectorCellIntegral.VALUE_VALUE);
		return List.of(mass);
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
		return x -> new CoordinateVector(x.getLength());
	}
	
	@Override
	protected Function<CoordinateVector, CoordinateVector> velocityBoundaryValues()
	{
		return x -> CoordinateVector.getUnitVector(2, 0)
		                            .mul(getTime());
	}
	
	public static void main(final String[] args)
	{
		final TaylorHoodSpace backGround = new TaylorHoodSpace(CoordinateVector.fromValues(0, 0),
		                                                       CoordinateVector.fromValues(1, 1),
		                                                       new IntCoordinates(10, 10));
		final CircleVectorSpace particle1 = new CircleVectorSpace(CoordinateVector.fromValues(0.2, 0.2),
		                                                          0.1,
		                                                          1);
		final CircleVectorSpace particle2 = new CircleVectorSpace(CoordinateVector.fromValues(0.5, 0.5),
		                                                          0.1,
		                                                          1);
		final CircleVectorSpace particle3 = new CircleVectorSpace(CoordinateVector.fromValues(0.2, 0.7),
		                                                          0.1,
		                                                          2);
		backGround.assembleCells();
		backGround.assembleFunctions(1);
		particle1.assembleCells();
		particle1.assembleFunctions(1);
		particle2.assembleCells();
		particle2.assembleFunctions(2);
		particle3.assembleCells();
		particle3.assembleFunctions(1);
		final DLM2Particles dlmElast2 = new DLM2Particles(0.02, 3, backGround, List.of(particle1, particle2,
		                                                                               particle3));
		
		dlmElast2.loop();
		final MixedPlot2DTime UpPlot = new MixedPlot2DTime(dlmElast2.pressureValues,
		                                                   dlmElast2.velocityValues,
		                                                   30,
		                                                   "Up1");
		UpPlot.addOverlay(new DistortedOverlay(dlmElast2.plotPoints,
		                                       particle1,
		                                       dlmElast2.getIterateHistory()
		                                                .slice(new IntCoordinates(0,
		                                                                          backGround.getShapeFunctions()
		                                                                                    .size()),
		                                                       new IntCoordinates(dlmElast2.getIterateHistory()
		                                                                                   .getRows(),
		                                                                          backGround.getShapeFunctions()
		                                                                                    .size() + particle1.getShapeFunctions()
		                                                                                                       .size())),
		                                       5));
		PlotWindow.addPlot(UpPlot);
		final MixedPlot2DTime UpPlot1 = new MixedPlot2DTime(dlmElast2.pressureValues,
		                                                    dlmElast2.velocityValues,
		                                                    30,
		                                                    "Up2");
		UpPlot1.addOverlay(new DistortedOverlay(dlmElast2.plotPoints,
		                                        particle2,
		                                        dlmElast2.getIterateHistory()
		                                                 .slice(new IntCoordinates(0,
		                                                                           backGround.getShapeFunctions()
		                                                                                     .size() + 2 * particle1.getShapeFunctions()
		                                                                                                            .size()),
		                                                        new IntCoordinates(dlmElast2.getIterateHistory()
		                                                                                    .getRows(),
		                                                                           backGround.getShapeFunctions()
		                                                                                     .size() + 2 * particle1.getShapeFunctions()
		                                                                                                            .size()
			                                                                           + particle2.getShapeFunctions()
			                                                                                      .size())),
		                                        5));
		PlotWindow.addPlot(UpPlot1);
		final MixedPlot2DTime UpPlot2 = new MixedPlot2DTime(dlmElast2.pressureValues,
		                                                    dlmElast2.velocityValues,
		                                                    30,
		                                                    "Up3");
		UpPlot2.addOverlay(new DistortedOverlay(dlmElast2.plotPoints,
		                                        particle3,
		                                        dlmElast2.getIterateHistory()
		                                                 .slice(new IntCoordinates(0,
		                                                                           backGround.getShapeFunctions()
		                                                                                     .size() + 2 * particle1.getShapeFunctions()
		                                                                                                            .size()
			                                                                           + 2 * particle2.getShapeFunctions()
			                                                                                          .size()),
		                                                        new IntCoordinates(dlmElast2.getIterateHistory()
		                                                                                    .getRows(),
		                                                                           backGround.getShapeFunctions()
		                                                                                     .size() + 2 * particle1.getShapeFunctions()
		                                                                                                            .size()
			                                                                           + 2 * particle2.getShapeFunctions()
			                                                                                          .size() + particle3.getShapeFunctions()
			                                                                                                             .size())),
		                                        5));
		PlotWindow.addPlot(UpPlot2);
	}
}
