import basic.PlotWindow;
import basic.ScalarFunction;
import basic.ScalarPlot2D;
import io.vavr.Tuple2;
import linalg.*;
import mixed.*;
import multigrid.MGPreconditionerSpace;
import multigrid.Smoother;
import schwarz.ColoredCartesianSchwarz;
import schwarz.DirectSolver;
import schwarz.SchwarzSmoother;
import tensorproduct.ContinuousTPShapeFunction;
import tensorproduct.ContinuousTPVectorFunction;
import tensorproduct.TPVectorCellIntegral;
import tensorproduct.TPVectorRightHandSideIntegral;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

import java.util.ArrayList;
import java.util.List;

public class StokesMGSchwarz
{
	public static void main(final String[] args)
	{
		final int refinements = 4;
		final TPVectorCellIntegral<ContinuousTPVectorFunction> gradGrad =
			new TPVectorCellIntegral<>(ScalarFunction.constantFunction(StokesReferenceSolution.reynolds),
			                           TPVectorCellIntegral.GRAD_GRAD);
		final MixedCellIntegral<TPCell, ContinuousTPShapeFunction, ContinuousTPVectorFunction, QkQkFunction>
			divValue =
			new MixedTPCellIntegral<>(ScalarFunction.constantFunction(-1),
			                          MixedTPCellIntegral.DIV_VALUE);
		final MixedCellIntegral<TPCell, ContinuousTPShapeFunction, ContinuousTPVectorFunction,
			QkQkFunction> vv
			= MixedCellIntegral.fromVelocityIntegral(gradGrad);
		final MixedRightHandSideIntegral<TPCell, ContinuousTPShapeFunction,
			ContinuousTPVectorFunction, QkQkFunction> rightHandSideIntegral =
			MixedRightHandSideIntegral.fromVelocityIntegral(
				new TPVectorRightHandSideIntegral<>(StokesReferenceSolution.rightHandSide(),
				                                    TPVectorRightHandSideIntegral.VALUE));
		final MGPreconditionerSpace<TaylorHoodSpace, TPCell, TPFace, QkQkFunction,
			MixedValue, MixedGradient, MixedHessian>
			mg = new MGPreconditionerSpace<>(refinements,
			                                 1)
		{

//			@Override
//			public void presmoothcallback(final int level, final Vector guess, final Vector rhs)
//			{
//				if (level == refinements)
//				{
//					final MixedTPFESpaceFunction<QkQkFunction> guessFUnction =
//						new MixedTPFESpaceFunction<>(getSpace(level).getShapeFunctionMap(),
//						                             guess);
//					PlotWindow.addPlot(new MixedPlot2D(guessFUnction,
//					                                   getFinestSpace().generatePlotPoints(40),
//					                                   40, "presmoot"));
//				}
//			}
//
//			@Override
//			public void precorrectioncallback(final int level, final Vector guess, final Vector rhs)
//			{
//				if (level == refinements)
//				{
//					final MixedTPFESpaceFunction<QkQkFunction> guessFUnction =
//						new MixedTPFESpaceFunction<>(getSpace(level).getShapeFunctionMap(),
//						                             guess);
//					PlotWindow.addPlot(new MixedPlot2D(guessFUnction,
//					                                   getFinestSpace().generatePlotPoints(40),
//					                                   40, "precorrect"));
//				}
//			}
//
//			@Override
//			public void postcorrectioncallback(final int level, final Vector guess,
//			                                   final Vector rhs)
//			{
//				if (level == refinements)
//				{
//					final MixedTPFESpaceFunction<QkQkFunction> guessFUnction =
//						new MixedTPFESpaceFunction<>(getSpace(level).getShapeFunctionMap(),
//						                             guess);
//					PlotWindow.addPlot(new MixedPlot2D(guessFUnction,
//					                                   getFinestSpace().generatePlotPoints(
//						                                   40),
//					                                   40, "postcorrect"));
//				}
//			}
			
			@Override
			public List<TaylorHoodSpace> createSpaces(final int refinements)
			{
				final ArrayList<TaylorHoodSpace> ret = new ArrayList<>();
				int mul = 1;
				for (int i = 0; i < refinements + 1; i++)
				{
					final TaylorHoodSpace s
						= new TaylorHoodSpace(CoordinateVector.fromValues(0, 0),
						                      CoordinateVector.fromValues(1, 1),
						                      new IntCoordinates(4,
						                                         4).mul(mul));
					ret.add(s);
					mul *= 2;
				}
				return ret;
			}
			
			@Override
			public Tuple2<VectorMultiplyable, DenseVector> createSystem(
				final TaylorHoodSpace space)
			{
				final SparseMatrix s = new SparseMatrix(space.getShapeFunctionMap()
				                                             .size(),
				                                        space.getShapeFunctionMap()
				                                             .size());
				final DenseVector rhs = new DenseVector(space.getShapeFunctionMap()
				                                             .size());
				space.writeCellIntegralsToMatrix(List.of(vv, divValue), s);
				space.writeCellIntegralsToRhs(List.of(rightHandSideIntegral), rhs);
				space.writeBoundaryValuesTo(new ComposedMixedFunction(StokesReferenceSolution.vectorBoundaryValues()),
				                            f -> true,
				                            (f, sf) -> sf.hasVelocityFunction(),
				                            s,
				                            rhs);
				final int firstPressure =
					space.getVelocitySize();
				space.overWriteValue(firstPressure, 0, s, rhs);
				return new Tuple2<>(s, rhs);
			}
			
			@Override
			public List<Smoother> createSmoothers()
			{
				final ArrayList<Smoother> ret = new ArrayList<>();
				for (int i = 1; i < spaces.size(); i++)
				{
					final IntCoordinates partitions
						= new IntCoordinates(2, 2).mul(Math.max(1, (int) Math.pow(2, i - 1)));
					System.out.println(partitions);
//					final VankaSchwarz schwarz =
//						new VankaSchwarz((SparseMatrix) getSystem(i),
//						                 getSpace(i),
//						                 new MultiplicativeSubspaceCorrection<>(getSpace(i)),
//						                 new DirectSolver());

//					final CartesianUpFrontSchwarz<QkQkFunction> schwarz
//						= new CartesianUpFrontSchwarz<>((SparseMatrix) getSystem(i),
//						                                getSpace(i),
//						                                partitions, 2,
//						                                new MultiplicativeSubspaceCorrection<>(),
//						                                new DirectSolver());
					final ColoredCartesianSchwarz<QkQkFunction> schwarz
						= new ColoredCartesianSchwarz<>((SparseMatrix) getSystem(i),
						                                getSpace(i),
						                                partitions, 2,
						                                new DirectSolver(), 1);
					ret.add(new SchwarzSmoother(2, schwarz));
				}
				return ret;
			}
			
			@Override
			public void applyZeroBoundaryConditions(final TaylorHoodSpace space,
			                                        final MutableVector vector)
			{
				space.projectOntoBoundaryValues(
					new ComposedMixedFunction(ScalarFunction.constantFunction(0)
					                                        .makeIsotropicVectorFunction()),
					f -> true,
					(f, sf) -> sf.hasVelocityFunction(),
					vector);
				vector.set(0, space.getVelocitySize());
			}
		};
		mg.verbose = true;
		
		PlotWindow.addPlot(new ScalarPlot2D(LaplaceReferenceSolution.scalarReferenceSolution(),
		                                    mg.spaces.get(refinements)
		                                             .generatePlotPoints(30),
		                                    30, "reference"));
		//		final IterativeSolver it = new IterativeSolver(true);
//		it.showProgress = true;
//		final Vector solution = it.solvePGMRES(mg.finest_system,
//		                                       mg,
//		                                       mg.finest_rhs,
//		                                       1e-8);
//		final VankaSchwarz schwarz =
//			new VankaSchwarz((SparseMatrix) mg.getFinestSystem(),
//			                 mg.getFinestSpace(),
//			                 new MultiplicativeSubspaceCorrection<>(mg.getFinestSpace()),
//			                 new DirectSolver());
//		final ColoredCartesianSchwarz<QkQkFunction> schwarz
//			= new ColoredCartesianSchwarz<>((SparseMatrix) mg.getFinestSystem(),
//			                                mg.getFinestSpace(),
//			                                new IntCoordinates(8, 8), 2,
//			                                new DirectSolver());
//		final CartesianUpFrontSchwarz<QkQkFunction> schwarz
//			= new CartesianUpFrontSchwarz<>((SparseMatrix) mg.getFinestSystem(),
//			                                mg.getFinestSpace(),
//			                                new IntCoordinates(8, 8), 2,
//			                                new MultiplicativeSubspaceCorrection<>(),
//			                                new DirectSolver());
		final Vector iterate = new DenseVector(mg.finest_rhs.getLength());
		final double initialres = mg.getFinestSystem()
		                            .mvMul(iterate)
		                            .sub(mg.finest_rhs)
		                            .euclidianNorm();
		final IterativeSolverConvergenceMetric m = new IterativeSolverConvergenceMetric(1e-8 * initialres);
		final IterativeSolver it = new IterativeSolver("it");
		it.showProgress = true;
		for (int i = 0; i < iterate.getLength(); i++)
			((DenseVector) iterate).
				set(Math.random() * 8 - 4, i);
		mg.applyZeroBoundaryConditions(mg.getFinestSpace(), ((DenseVector) iterate));
		System.out.println("SOLVING SYSTEM OF SIZE " + mg.getFinestSystem()
		                                                 .getTVectorSize());
//		MetricWindow.getInstance()
//		            .setMetric("it", m);
//		for (int i = 0; i < 100; i++)
//		{
//			iterate = schwarz.getSubspaceCorrection()
//			                 .apply(schwarz, iterate, mg.finest_rhs);
//			m.publishIterate(mg.getFinestSystem()
//			                   .mvMul(iterate)
//			                   .sub(mg.finest_rhs)
//			                   .euclidianNorm());
//			System.out.println(mg.getFinestSystem()
//			                     .mvMul(iterate)
//			                     .sub(mg.finest_rhs)
//			                     .euclidianNorm());
//		}
//		iterate = it.solvePGMRES(mg.finest_system, schwarz, mg.finest_rhs, 1e-6);
//		double res = initialres;
//		MetricWindow.getInstance()
//		            .setMetric("it", m);
//		for (int i = 0; i < 100 && res > initialres * 1e-8; i++)
//		{
//			iterate = mg.vCycle(iterate, mg.finest_rhs);
//			res = mg.finest_system.mvMul(iterate)
//			                      .sub(mg.finest_rhs)
//			                      .euclidianNorm();
//			m.publishIterate(res);
//			System.out.println(res / initialres + " after iterations " + i);
//		}
		it.solvePGMRES(mg.finest_system, mg, mg.finest_rhs, 1e-8);
		System.out.println("SOLVED");
	}
}
