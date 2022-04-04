import basic.*;
import io.vavr.Tuple2;
import linalg.Vector;
import linalg.*;
import mixed.*;
import multigrid.MGPreconditionerSpace;
import multigrid.Smoother;
import org.jetbrains.annotations.NotNull;
import schwarz.CartesianUpFrontSchwarz;
import schwarz.DirectSolver;
import schwarz.MultiplicativeSubspaceCorrection;
import schwarz.SchwarzSmoother;
import tensorproduct.*;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

import java.awt.*;
import java.util.List;
import java.util.*;

public class LidDriven
{
	
	private static MGPreconditionerSpace<TaylorHoodSpace, TPCell, TPFace, QkQkFunction, MixedValue, MixedGradient
		, MixedHessian> getMG(final Vector guess)
	{
		final double reyn = 100;
		final MixedCellIntegral<TPCell, ContinuousTPShapeFunction, ContinuousTPVectorFunction, QkQkFunction>
			divValue =
			new MixedTPCellIntegral<>(ScalarFunction.constantFunction(1),
			                          MixedTPCellIntegral.DIV_VALUE);
		final MixedCellIntegral<TPCell, ContinuousTPShapeFunction, ContinuousTPVectorFunction, QkQkFunction>
			gradGrad =
			MixedCellIntegral.fromVelocityIntegral(new TPVectorCellIntegral<>(
				ScalarFunction.constantFunction(1. / reyn),
				TPVectorCellIntegral.GRAD_GRAD));
		return new MGPreconditionerSpace<>(1, 1)
		{
			
			@Override
			public List<TaylorHoodSpace> createSpaces(final int refinements)
			{
				verbose = true;
				final ArrayList<TaylorHoodSpace> ret = new ArrayList<>();
				int mul = 1;
				for (int i = 0; i < refinements + 1; i++)
				{
					final TaylorHoodSpace s = new TaylorHoodSpace(CoordinateVector.fromValues(0,
					                                                                          0),
					                                              CoordinateVector.fromValues(1,
					                                                                          1),
					                                              new IntCoordinates(8,
					                                                                 8).mul(mul));
					ret.add(s);
					mul *= 2;
				}
				return ret;
			}
			
			Map<Integer, Integer> sizeToFixedP;
			
			@Override
			public Tuple2<VectorMultiplyable, DenseVector> createSystem(final TaylorHoodSpace space)
			{
				if (sizeToFixedP == null)
					sizeToFixedP = new HashMap<>();
				final int n = space.getShapeFunctionMap()
				                   .size();
				final SparseMatrix flow = new SparseMatrix(n, n);
				final List<CellIntegral<TPCell, QkQkFunction>> flowIntegrals = new ArrayList<>();
				if (guess != null)
				{
					flowIntegrals.addAll(
						getConvectionIntegrals(generateCurrentFunction(guess, getFinestSpace())
							                       .getVelocityFunction()));
				}
				flowIntegrals.add(gradGrad);
				flowIntegrals.add(divValue);
				space.writeCellIntegralsToMatrix(flowIntegrals, flow);
				final DenseVector src = new DenseVector(n);
				space.writeCellIntegralsToRhs(List.of(getSourceIntegral()), src);
				
				space.writeBoundaryValuesTo(getBoundaryFunction(),
				                            f -> true,
				                            (f, sf) -> sf.hasVelocityFunction(),
				                            flow,
				                            src);
				final int fixedP = getFixedP(space);
				
				final DenseVector cons = new DenseVector(n);
				for (int i = space.getVelocitySize(); i < flow.getCols(); i++)
					cons.set(1, i);
//				final MixedTPFESpaceFunction<QkQkFunction> constantFunction
//					= new MixedTPFESpaceFunction<>(
//					space.getShapeFunctionMap(),
//					cons);
//				PlotWindow.addPlot(new MixedPlot2D(constantFunction,
//				                                   space.generatePlotPoints(30),
//				                                   30));
				System.out.println("COOOOOOOOOOOOOOOONNSSSS " + flow.mvMul(cons)
				                                                    .euclidianNorm());
				space.overWriteValue(fixedP, 0, flow, src);
				sizeToFixedP.put(n, fixedP);
				return new Tuple2<>(flow, src);
			}
			
			private int getFixedP(final TaylorHoodSpace space)
			{
				int fixedP = 0;
				final CoordinateVector c =
					space.grid.startCoordinates.add(space.grid.endCoordinates)
					                           .mul(0.5);
				fixedP = space.getShapeFunctionMap()
				              .entrySet()
				              .stream()
				              .filter(e -> e.getValue()
				                            .hasPressureFunction())
				              .min(Comparator.comparing(e -> e.getValue()
				                                              .getPressureFunction()
				                                              .getNodeFunctional()
				                                              .getPoint()
				                                              .dist(c)))
				              .orElseThrow()
				              .getKey();
				return space.getVelocitySize();//fixedP;
			}
			
			@Override
			public List<Smoother> createSmoothers()
			{
				final ArrayList<Smoother> ret = new ArrayList<>();
				for (int i = 1; i < spaces.size(); i++)
				{
					final IntCoordinates partitions
						= new IntCoordinates(1, 1).mul(Math.max(1, (int) Math.pow(2, i)));
					System.out.println(partitions);
					final CartesianUpFrontSchwarz<QkQkFunction> schwarz =
						new CartesianUpFrontSchwarz<>((SparseMatrix) getSystem(i),
						                              getSpace(i),
						                              partitions,
						                              1,
						                              new MultiplicativeSubspaceCorrection<>(
							                              getSpace(i)),
						                              new DirectSolver());
					ret.add(new SchwarzSmoother(3, schwarz));
				}
				return ret;
			}
			
			@Override
			public void applyZeroBoundaryConditions(final TaylorHoodSpace space, final MutableVector vector)
			{
				space.projectOntoBoundaryValues(new ComposedMixedFunction(ScalarFunction.constantFunction(
					                                                                        0)
				                                                                        .makeIsotropicVectorFunction()),
				                                f -> true,
				                                (f, sf) -> sf.hasVelocityFunction(),
				                                vector);
				if (sizeToFixedP.containsKey(vector.getLength()))
					vector.set(0, sizeToFixedP.get(vector.getLength()));
			}
		};
	}
	
	private static MixedFunction getBoundaryFunction()
	{
		return new ComposedMixedFunction(new VectorFunction()
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
				if (pos.y() == 1)
					return CoordinateVector.fromValues(1, 0);
				return new CoordinateVector(2);
			}
		});
	}
	
	public static MixedFunctionOnCells<TPCell, TPFace> generateCurrentFunction(final Vector iterate,
	                                                                           final TaylorHoodSpace grid)
	{
		return new MixedTPFESpaceFunction<>(grid.getShapeFunctionMap(), iterate);
	}
	
	public static void main(final String[] args)
	{
		
		final PerformanceArguments.PerformanceArgumentBuilder builder =
			new PerformanceArguments.PerformanceArgumentBuilder();
		builder.build();
		final int nPoints = 80;
		DenseVector guess = null;
		final IterativeSolverConvergenceMetric cm = new IterativeSolverConvergenceMetric(1e-10);
		MetricWindow.getInstance()
		            .setMetric("nonlinear", cm);
		MGPreconditionerSpace<TaylorHoodSpace, TPCell, TPFace, QkQkFunction, MixedValue, MixedGradient,
			MixedHessian> mg = null;
		for (int nLinIter = 0; nLinIter < 20; nLinIter++)
		{
			final IterativeSolverConvergenceMetric cm2 = new IterativeSolverConvergenceMetric(1e-10);
			MetricWindow.getInstance()
			            .setMetric("itsolver", cm2);
			mg = getMG(guess);
			final DenseMatrix d = new DenseMatrix((SparseMatrix) mg.finest_system);
			System.out.println("copied" + d.getShape());
			
			System.out.println("SIZE " + mg.finest_system.getVectorSize());
			if (guess == null)
				guess = new DenseVector(mg.fullVCycleCorrection(mg.finest_rhs));
			mg.applyBoundaryConditions(mg.getFinestSpace(), guess, mg.finest_rhs);
			final Vector b = mg.finest_rhs;
			final Vector defect = b.sub(mg.getFinestSystem()
			                              .mvMul(guess));
			System.out.println("NLIN DEFECT " + defect.euclidianNorm());
			Vector correct = new DenseVector(guess.getLength());
			final double initialres = mg.getFinestSystem()
			                            .mvMul(correct)
			                            .sub(defect)
			                            .euclidianNorm();
			double res = initialres;
			int iter;
			//correct = ((SparseMatrix) mg.finest_system).solveNative(defect);
			if (mg.finest_system.getVectorSize() < 2000)
				if (((SparseMatrix) mg.finest_system).inverseNative()
				                                     .absMaxElement() > 1e10)
					throw new IllegalStateException("System is not invertivle");
			System.out.println(defect.sub(mg.getFinestSystem()
			                                .mvMul(correct))
			                         .euclidianNorm());
			for (iter = 0; iter < 100 && res > 1e-8; iter++)
			
			{
				cm2.publishIterate(res);
				correct = new DenseVector(mg.vCycle(correct, defect));
				res = mg.finest_system.mvMul(correct)
				                      .sub(defect)
				                      .euclidianNorm();
				System.out.println("NLIN ITER " + nLinIter + " " + res / initialres + " after " +
					                   "iterations " + iter);
				System.out.println("lin norm " + mg.finest_system.mvMul(correct)
				                                                 .sub(defect)
				                                                 .euclidianNorm());
			}
			if (iter > 90)
			{
				break;
			}
			//correct = new IterativeSolver(true).solvePGMRES(mg.getFinestSystem(), mg, defect, 1e-8);
			cm.publishIterate(correct.euclidianNorm());
			System.out.println("GUESS ITERATE" + nLinIter);
			System.out.println(defect.euclidianNorm());
			System.out.println("NONLINEARN CORRECTION " + correct.euclidianNorm() + " " + guess.euclidianNorm());
			guess = guess.add(correct);
			if (correct.euclidianNorm() < 1e-8)
				break;
		}
		final Plot p = new VelocityMagnitudePlot2D(generateCurrentFunction(guess,
		                                                                   mg.getFinestSpace()),
		                                           mg.getFinestSpace()
		                                             .generatePlotPoints(nPoints), nPoints
			, "FEM");
		final CoordinateSystemOverlay o
			= new CoordinateSystemOverlay(mg.getFinestSpace().grid.startCoordinates,
			                              mg.getFinestSpace().grid.endCoordinates);
		for (final var f : mg.getFinestSpace()
		                     .getShapeFunctions())
		{
			if (f.hasVelocityFunction())
			{
				o.addPoint(f.getVelocityFunction()
				            .getNodeFunctionalPoint(),
				           "V" + f.getVelocityFunction()
				                  .getComponent(), 10, Color.GREEN);
			}
			if (f.hasPressureFunction())
			{
				o.addPoint(f.getPressureFunction()
				            .getNodeFunctional()
				            .getPoint(), "P   ", 6, Color.RED);
			}
		}
		p.addOverlay(o);
		PlotWindow.addPlotShow(p);
		PlotWindow.addPlot(new MixedPlot2D(generateCurrentFunction(guess, mg.getFinestSpace()),
		                                   mg.getFinestSpace()
		                                     .generatePlotPoints(80)
			,
			                           80,
			                           "fem"));
	}
	
	@NotNull
	private static MixedRightHandSideIntegral<TPCell, ContinuousTPShapeFunction, ContinuousTPVectorFunction,
		QkQkFunction> getSourceIntegral()
	{
		return MixedRightHandSideIntegral.fromVelocityIntegral(
			new TPVectorRightHandSideIntegral<>(
				ScalarFunction.constantFunction(0)
				              .makeIsotropicVectorFunction(),
				TPVectorRightHandSideIntegral.VALUE));
	}
	
	public static List<CellIntegral<TPCell, QkQkFunction>> getConvectionIntegrals(final VectorFunctionOnCells<TPCell,
		TPFace> velocity)
	{
		final VectorFunctionOnCells<TPCell, TPFace> semiImplicitWeight1 =
			VectorFunctionOnCells.fromLambda(velocity::value,
			                                 velocity::valueInCell, 2, 2);
		final TPVectorCellIntegralOnCell<ContinuousTPVectorFunction> convection1 =
			new TPVectorCellIntegralOnCell<>(semiImplicitWeight1, TPVectorCellIntegral.GRAD_VALUE);
		final MixedCellIntegral<TPCell, ContinuousTPShapeFunction, ContinuousTPVectorFunction,
			QkQkFunction> mixedConvection1 = MixedCellIntegral.fromVelocityIntegral(convection1);
		return List.of(mixedConvection1);
	}
}
