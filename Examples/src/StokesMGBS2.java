import basic.PlotWindow;
import basic.ScalarFunction;
import basic.ScalarPlot2D;
import dlm.BSSmoother5;
import io.vavr.Tuple2;
import linalg.*;
import mixed.*;
import multigrid.MGPreconditionerSpace;
import multigrid.Smoother;
import tensorproduct.ContinuousTPShapeFunction;
import tensorproduct.ContinuousTPVectorFunction;
import tensorproduct.TPVectorCellIntegral;
import tensorproduct.TPVectorRightHandSideIntegral;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

import java.util.ArrayList;
import java.util.List;

public class StokesMGBS2
{
	public static void main(final String[] args)
	{
		final int refinements = 3;
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
					ret.add(new BSSmoother5(3, 6, 1, getSpace(i).getVelocitySize(), this));
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
		Vector iterate = new DenseVector(mg.finest_rhs.getLength());
		for (int i = 0; i < iterate.getLength(); i++)
			((DenseVector) iterate).
				
				set(Math.random(), i);
		mg.applyZeroBoundaryConditions(mg.getFinestSpace(), ((DenseVector) iterate));
		final double initialres = mg.getFinestSystem()
		                            .mvMul(iterate)
		                            .sub(mg.finest_rhs)
		                            .euclidianNorm();
		double res = initialres;
		for (int i = 0; i < 100 && res > initialres * 1e-8; i++)
		
		{
			iterate = mg.vCycle(iterate, mg.finest_rhs);
			res = mg.finest_system.mvMul(iterate)
			                      .sub(mg.finest_rhs)
			                      .euclidianNorm();
			System.out.println(res / initialres + " after iterations " + i);
		}
	}
}
