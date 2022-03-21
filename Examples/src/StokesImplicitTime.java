import basic.*;
import io.vavr.Tuple2;
import linalg.*;
import mixed.*;
import multigrid.MGPreconditionerSpace;
import multigrid.Smoother;
import org.jetbrains.annotations.NotNull;
import schwarz.ColoredCartesianSchwarz;
import schwarz.DirectSolver;
import schwarz.SchwarzSmoother;
import tensorproduct.*;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

import java.util.ArrayList;
import java.util.List;

public class StokesImplicitTime
{
	public static MGPreconditionerSpace<TaylorHoodSpace, TPCell, TPFace, QkQkFunction, MixedValue, MixedGradient
		, MixedHessian> getMG(final double dt, final double t, final Vector previous, final Vector guess)
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
		return new MGPreconditionerSpace<>(2, 1)
		{
			
			@Override
			public List<TaylorHoodSpace> createSpaces(final int refinements)
			{
				verbose = false;
				final ArrayList<TaylorHoodSpace> ret = new ArrayList<>();
				int mul = 1;
				for (int i = 0; i < refinements + 1; i++)
				{
					final TaylorHoodSpace s = new TaylorHoodSpace(CoordinateVector.fromValues(0, 0),
					                                              CoordinateVector.fromValues(1, 1),
					                                              new IntCoordinates(4,
					                                                                 4).mul(mul));
					ret.add(s);
					mul *= 2;
				}
				return ret;
			}
			
			@Override
			public Tuple2<VectorMultiplyable, DenseVector> createSystem(final TaylorHoodSpace space)
			{
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
						getConvectionIntegrals(generateCurrentFunction(guess,
						                                               getFinestSpace())
							                       .getVelocityFunction()));
				}
				flowIntegrals.add(gradGrad);
				flowIntegrals.add(divValue);
				space.writeCellIntegralsToMatrix(flowIntegrals, flow);
				final DenseVector src = new DenseVector(n);
				space.writeCellIntegralsToRhs(List.of(getSourceIntegral(t)), src);
				
				final SparseMatrix s = mass.add(flow);
				DenseVector fromLast = new DenseVector(n);
				if (previous != null && previous.getLength() == n)
					fromLast = mass.mvMul(previous);
				final DenseVector d = src.add(fromLast);
				
				space.writeBoundaryValuesTo(new ComposedMixedFunction(getBoundaryFunction(t)),
				                            f -> true,
				                            (f, sf) -> sf.hasVelocityFunction(),
				                            s,
				                            d);
				
				space.overWriteValue(space.getVelocitySize(), 0, s, d);
				return new Tuple2<>(s, d);
			}
			
			@Override
			public List<Smoother> createSmoothers()
			{
				final ArrayList<Smoother> ret = new ArrayList<>();
				for (int i = 1; i < spaces.size(); i++)
				{
					final IntCoordinates partitions
						= new IntCoordinates(2, 2).mul(Math.max(1, (int) Math.pow(2, i)));
					System.out.println(partitions);
//					final VankaSchwarz schwarz =
//						new VankaSchwarz((SparseMatrix) getSystem(i),
//						                 getSpace(i),
//						                 new MultiplicativeSubspaceCorrection<>(getSpace(i)),
//						                 new DirectSolver());
					
					final ColoredCartesianSchwarz<QkQkFunction> schwarz
						= new ColoredCartesianSchwarz<>(
						(SparseMatrix) getSystem(i),
						getSpace(i),
						partitions, 1,
						new DirectSolver(), 0.1);
					ret.add(new SchwarzSmoother(2, schwarz));
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
				vector.set(0, space.getVelocitySize());
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
				final double effectiveT = Math.min(1, t);
				if (Math.abs(pos.x()) <= 1e-10 || Math.abs(pos.x() - 1) <= 1e-10)
					return CoordinateVector.fromValues(effectiveT * .3, effectiveT * 1.5);
				///if(Math.abs(pos.x()-1) <= 1e-1)
				//	return CoordinateVector.fromValues(1,1).mul(Math.sin(t)*Math.sin(t)*10*(0.25-
				//	(0.5-pos.y())*(0.5-pos.y())));
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
		final int nPoints = 60;
		
		MGPreconditionerSpace<TaylorHoodSpace, TPCell, TPFace, QkQkFunction, MixedValue, MixedGradient
			, MixedHessian> mg = getMG(dt, 0, null, null);
		
		Vector iterate = getInitialIterate(mg.getFinestSpace());
		
		final List<CoordinateVector> points = mg.getFinestSpace()
		                                        .generatePlotPoints(nPoints);
		
		for (int i = 1; i <= timesteps; i++)
		{
			Vector guess = new DenseVector(iterate);
//			final IterativeSolverConvergenceMetric nm =
//				new IterativeSolverConvergenceMetric(guess.euclidianNorm() * 1e-5);
//			MetricWindow.getInstance()
//			            .setMetric("nlin", nm);
			for (int nLinIter = 0; nLinIter < 10; nLinIter++)
			{
				mg = getMG(dt, i * dt, iterate, guess);
				final Vector b = mg.finest_rhs;
				final Vector defect = b.sub(mg.getFinestSystem()
				                              .mvMul(guess));
				System.out.println("NLIN DEFECT " + defect.euclidianNorm());
				//final IterativeSolverConvergenceMetric cm = new IterativeSolverConvergenceMetric
				// (1e-8);
//				MetricWindow.getInstance()
//				            .setMetric("mg", cm);
				final Vector correct = new IterativeSolver("mg").solvePGMRES(mg.getFinestSystem(),
				                                                             mg,
				                                                             defect,
				                                                             1e-8);
//				Vector correct = guess.mul(0);
//				final double initialres = mg.getFinestSystem()
//				                            .mvMul(correct)
//				                            .sub(defect)
//				                            .euclidianNorm();
//				double res = initialres;
//				int it;
//				for (it = 0; it < 100 && res > 1e-8; it++)
//
//				{
//					cm.publishIterate(res);
//					correct = mg.vCycle(correct, defect);
//					res = mg.finest_system.mvMul(correct)
//					                      .sub(defect)
//					                      .euclidianNorm();
//					System.out.println("NLIN ITER " + nLinIter + " " + res / initialres + " after " +
//						                   "iterations " + it);
//				}
//				if (it > 30)
//				{
//					i = 2 * timesteps + 1;
//					break;
//				}
//				nm.publishIterate(correct.euclidianNorm());
				System.out.println("NONLINEARN CORRECTION " + correct.euclidianNorm() + " " + guess.euclidianNorm() + " " + iterate.euclidianNorm());
				guess = guess.add(correct.mul(1));
				if (correct.euclidianNorm() < 1e-5 * iterate.euclidianNorm())
					break;
			}
			iterate = guess;
			PlotWindow.addPlotShow(new MixedPlot2D(generateCurrentFunction(iterate,
			                                                               mg.getFinestSpace()),
			                                       points, nPoints, "TIME " + i * dt));
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
						return CoordinateVector.fromValues(0,
						                                   0 * 1 * Math.sin((pos.x() - time) * 10));
					}
				},
				TPVectorRightHandSideIntegral.VALUE));
	}
	
	public static List<CellIntegral<TPCell, QkQkFunction>> getConvectionIntegrals(final VectorFunctionOnCells<TPCell,
		TPFace> velocity)
	{
		
		final VectorFunctionOnCells<TPCell, TPFace> semiImplicitWeight1 =
			VectorFunctionOnCells.fromLambda((x) -> velocity.value(x)
			                                                .mul(-1. / 2),
			                                 (x, cell) ->
				                                 velocity.valueInCell(x, cell)
				                                         .mul(1. / 2), 2, 2);
		final VectorFunctionOnCells<TPCell, TPFace> semiImplicitWeight2 =
			VectorFunctionOnCells.fromLambda((x) -> velocity.value(x)
			                                                .mul(1. / 2),
			                                 (x, cell) -> velocity.valueInCell(x, cell)
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
}
