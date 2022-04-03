import basic.*;
import io.vavr.Tuple2;
import linalg.Vector;
import linalg.*;
import mixed.*;
import multigrid.MGPreconditionerSpace;
import multigrid.Smoother;
import org.jetbrains.annotations.NotNull;
import schwarz.DirectSolver;
import schwarz.MultiplicativeSubspaceCorrection;
import schwarz.SchwarzSmoother;
import schwarz.VankaSchwarz;
import tensorproduct.*;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

import java.awt.*;
import java.util.List;
import java.util.*;

public class StokesMGVankaTime
{
	private static MGPreconditionerSpace<TaylorHoodSpace, TPCell, TPFace, QkQkFunction, MixedValue, MixedGradient
		, MixedHessian> getMG(final Vector previous, final Vector guess, final double dt, final double time)
	{
		final double reynolds = 1000;
		final MixedCellIntegral<TPCell, ContinuousTPShapeFunction, ContinuousTPVectorFunction, QkQkFunction>
			divValue =
			new MixedTPCellIntegral<>(ScalarFunction.constantFunction(1),
			                          MixedTPCellIntegral.DIV_VALUE);
		final MixedCellIntegral<TPCell, ContinuousTPShapeFunction, ContinuousTPVectorFunction, QkQkFunction>
			gradGrad =
			MixedCellIntegral.fromVelocityIntegral(new TPVectorCellIntegral<>(
				ScalarFunction.constantFunction(1. / reynolds),
				TPVectorCellIntegral.GRAD_GRAD));
		final MixedCellIntegral<TPCell, ContinuousTPShapeFunction, ContinuousTPVectorFunction, QkQkFunction>
			valueValue =
			MixedCellIntegral.fromVelocityIntegral(new TPVectorCellIntegral<>(
				ScalarFunction.constantFunction(1),
				TPVectorCellIntegral.VALUE_VALUE));
		return new MGPreconditionerSpace<>(1, 1)
		{

//			@Override
//			public void postmoothcallback(final int level, final Vector guess, final Vector rhs)
//			{
//				if (level == 0)
//				{
//					final var function =
//						new MixedTPFESpaceFunction<>(getSpace(level).getShapeFunctionMap(),
//						                             guess);
//					PlotWindow.addPlot(new MixedPlot2D(function,
//					                                   getSpace(level).generatePlotPoints(40), 40
//						, "coarse correction"));
//				}
//			}
//
//			@Override
//			public void precorrectioncallback(final int level, final Vector guess, final Vector rhs)
//			{
//				if (level == 1)
//				{
//					final Vector realCorrect =
//						((SparseMatrix) getSystem(level)).solve(
//							rhs.sub(getSystem(level).mvMul(guess)));
//					final var function =
//						new MixedTPFESpaceFunction<>(getSpace(level).getShapeFunctionMap(),
//						                             realCorrect);
//					PlotWindow.addPlot(new MixedPlot2D(function,
//					                                   getSpace(level).generatePlotPoints(40), 40
//						, "real pre correction"));
//				}
//			}
//
//			@Override
//			public void correctioncallback(final int level, final Vector correction, final Vector rhs)
//			{
//				if (level == 1)
//				{
//					final var function =
//						new MixedTPFESpaceFunction<>(getSpace(level).getShapeFunctionMap(),
//						                             correction);
//					PlotWindow.addPlot(new MixedPlot2D(function,
//					                                   getSpace(level).generatePlotPoints(40), 40
//						, "actual correction"));
//				}
//			}
			
			@Override
			public List<TaylorHoodSpace> createSpaces(final int refinements)
			{
				verbose = true;
				final ArrayList<TaylorHoodSpace> ret = new ArrayList<>();
				int mul = 1;
				for (int i = 0; i < refinements + 1; i++)
				{
					final TaylorHoodSpace s = new TaylorHoodSpace(CoordinateVector.fromValues(0, 0),
					                                              CoordinateVector.fromValues(1, 1),
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
				final SparseMatrix mass = new SparseMatrix(n, n);
				space.writeCellIntegralsToMatrix(List.of(valueValue), mass);
				mass.mulInPlace(1. / dt);
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
				space.writeCellIntegralsToRhs(List.of(getSourceIntegral(time)), src);
				
				final SparseMatrix s = mass.add(flow);
				DenseVector fromLast = new DenseVector(n);
				if (previous != null && previous.getLength() == n)
					fromLast = mass.mvMul(previous);
				final DenseVector d = src.add(fromLast);
				
				space.writeBoundaryValuesTo(new ComposedMixedFunction(getBoundaryFunction(time)),
				                            f -> true,
				                            (f, sf) -> sf.hasVelocityFunction(),
				                            s,
				                            d);
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
				
				final DenseVector cons = new DenseVector(n);
				for (int i = space.getVelocitySize(); i < s.getCols(); i++)
					cons.set(1, i);
				final MixedTPFESpaceFunction<QkQkFunction> constantFunction
					= new MixedTPFESpaceFunction<>(
					space.getShapeFunctionMap(),
					cons);
//				PlotWindow.addPlot(new MixedPlot2D(constantFunction,
//				                                   space.generatePlotPoints(30),
//				                                   30));
				System.out.println("COOOOOOOOOOOOOOOONNSSSS " + s.mvMul(cons)
				                                                 .euclidianNorm());
				cons.set(fixedP, 0);
				System.out.println("COOOOOOOOOOOOOOOONNSSSS " + s.mvMul(cons)
				                                                 .euclidianNorm());
//				for (int i = space.getVelocitySize(); i < s.getCols(); i++)
//				{
//					s.set(0, i, fixedP);
//					s.set(0, fixedP, i);
//				}
				//d.set(0, fixedP);
				sizeToFixedP.put(n, fixedP);
				space.overWriteValue(fixedP, 0, s, d);
				return new Tuple2<>(s, d);
			}
			
			@Override
			public List<Smoother> createSmoothers()
			{
				final ArrayList<Smoother> ret = new ArrayList<>();
				for (int i = 1; i < spaces.size(); i++)
				{
					final IntCoordinates partitions
						= new IntCoordinates(4, 4).mul(Math.max(1, (int) Math.pow(2, i)));
					System.out.println(partitions);
					final VankaSchwarz schwarz =
						new VankaSchwarz((SparseMatrix) getSystem(i),
						                 getSpace(i),
						                 new MultiplicativeSubspaceCorrection<>(getSpace(i)),
						                 new DirectSolver());

//					final ColoredCartesianSchwarz<QkQkFunction> schwarz
//						= new ColoredCartesianSchwarz<>((SparseMatrix) getSystem(i),
//						                                getSpace(i),
//						                                partitions, 1,
//						                                new DirectSolver(), 1);
					ret.add(new SchwarzSmoother(4, schwarz));
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
				vector.set(0, sizeToFixedP.get(vector.getLength()));
				//vector.set(0, sizeToFixedP.get(vector.getLength()));//, space.getVelocitySize());
			}
		};
	}
	
	public static VectorFunction getBoundaryFunction(final double t)
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
				if (Math.abs(pos.x() * pos.y() * (1 - pos.x()) * (1 - pos.y())) > 1e-8)
					return CoordinateVector.fromValues(100);
				if (Math.abs(pos.y() - 1) < 1e-5)
					return CoordinateVector.getUnitVector(2, 0);
				return new CoordinateVector(2);
			}
		};
	}
	
	public static VectorFunction getInitialFunction()
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
	
	private static Vector getInitialIterate(final TaylorHoodSpace space)
	{
		final DenseVector initial = new DenseVector(space.getShapeFunctions()
		                                                 .size());
		final MixedFunction initialVelo
			= new ComposedMixedFunction(getInitialFunction());
		space.getShapeFunctionMap()
		     .values()
		     .forEach(function ->
		              {
			              initial.set(function.getNodeFunctional()
			                                  .evaluate(initialVelo), function.getGlobalIndex());
		              });
		return initial;
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
		final double dt = 0.1;
		final int timesteps = 400;
		final int nPoints = 61;
		
		MGPreconditionerSpace<TaylorHoodSpace, TPCell, TPFace, QkQkFunction, MixedValue, MixedGradient
			, MixedHessian> mg = getMG(null, null, dt, 0);
		
		Vector iterate = getInitialIterate(mg.getFinestSpace());
		
		final List<CoordinateVector> points = mg.getFinestSpace()
		                                        .generatePlotPoints(nPoints);
		//final IterativeSolver it = new IterativeSolver(true);
		for (int i = 1; i <= timesteps; i++)
		{
			DenseVector guess = new DenseVector(iterate);
//			final IterativeSolverConvergenceMetric nm =
//				new IterativeSolverConvergenceMetric(guess.euclidianNorm() * 1e-5);
//			MetricWindow.getInstance()
//			            .setMetric("nlin", nm);
			final IterativeSolverConvergenceMetric cm = new IterativeSolverConvergenceMetric(1e-10);
			MetricWindow.getInstance()
			            .setMetric("nonlinear", cm);
			for (int nLinIter = 0; nLinIter < 20; nLinIter++)
			{
				mg = getMG(iterate, guess, dt, i * dt);
				//mg.applyBoundaryConditions(mg.getFinestSpace(), guess, mg.finest_rhs);
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
				int it;
//				correct = ((SparseMatrix) mg.finest_system).solve(defect);
				for (it = 0; it < 100 && res > 1e-8; it++)
				
				{
					cm.publishIterate(res);
					correct = mg.vCycle(correct, defect);
					res = mg.finest_system.mvMul(correct)
					                      .sub(defect)
					                      .euclidianNorm();
					System.out.println("NLIN ITER " + nLinIter + " " + res / initialres + " after " +
						                   "iterations " + it);
					System.out.println("lin norm " + mg.finest_system.mvMul(correct)
					                                                 .sub(defect)
					                                                 .euclidianNorm());
				}
				if (it > 90)
				{
					i = 2 * timesteps + 1;
					break;
				}
//				correct = it.solvePGMRES(mg.getFinestSystem(), mg, defect, 1e-8);
				cm.publishIterate(correct.euclidianNorm());
				System.out.println("GUESS ITERATE" + nLinIter);
				System.out.println(defect.euclidianNorm());
				System.out.println("NONLINEARN CORRECTION " + correct.euclidianNorm() + " " + guess.euclidianNorm() + " " + iterate.euclidianNorm());
				guess = guess.add(correct);
				if (correct.euclidianNorm() < 1e-5 * iterate.euclidianNorm() || correct.euclidianNorm() < 1e-8)
					break;
			}
			iterate = guess;
			final Plot p = new VelocityMagnitudePlot2D(generateCurrentFunction(iterate,
			                                                                   mg.getFinestSpace()),
			                                           points, nPoints, "Time " + i * dt);
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
//			PlotWindow.addPlot(new MixedPlot2D(generateCurrentFunction(guess,
//			                                                           mg.getFinestSpace()),
//			                                   points, nPoints, "TIME " + i * dt));
			System.out.println("Time " + i * dt + " out of " + timesteps * dt + " done");
		}
	}
	
	@NotNull
	private static MixedRightHandSideIntegral<TPCell, ContinuousTPShapeFunction, ContinuousTPVectorFunction,
		QkQkFunction> getSourceIntegral(final double time)
	{
		return MixedRightHandSideIntegral.fromVelocityIntegral(
			new TPVectorRightHandSideIntegral<>(
				new VectorFunction()
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
						return CoordinateVector.fromValues(0, 0);
					}
				},
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
