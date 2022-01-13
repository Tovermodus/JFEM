import basic.PlotWindow;
import basic.ScalarFESpaceFunction;
import basic.ScalarFunction;
import basic.ScalarPlot2D;
import io.vavr.Tuple2;
import linalg.*;
import multigrid.AMGSpace;
import multigrid.MGPReconditioner;
import multigrid.RichardsonSmoother;
import multigrid.Smoother;
import tensorproduct.ContinuousTPFESpace;
import tensorproduct.ContinuousTPShapeFunction;
import tensorproduct.TPCellIntegral;
import tensorproduct.TPRightHandSideIntegral;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

import java.util.ArrayList;
import java.util.List;

public class LaplaceAMG
{
	public static void main(final String[] args)
	{
		final int refinements = 6;
		final AMGSpace<ContinuousTPFESpace, TPCell, TPFace, ContinuousTPShapeFunction, Double, CoordinateVector,
			CoordinateMatrix>
			mg = new AMGSpace<>(CoordinateVector.fromValues(0, 0),
			                    CoordinateVector.fromValues(1, 1),
			                    new IntCoordinates(4, 4),
			                    refinements,
			                    1)
		{
			@Override
			public List<ContinuousTPFESpace> createSpaces(final int refinements)
			{
				final ArrayList<ContinuousTPFESpace> ret = new ArrayList<>();
				int mul = 1;
				for (int i = 0; i < refinements + 1; i++)
				{
					ret.add(new ContinuousTPFESpace(startCoordinates, endCoordinates,
					                                coarseCellsPerDimension.mul(mul)));
					mul *= 2;
				}
				return ret;
			}
			
			@Override
			public Tuple2<SparseMatrix, Vector> createFinestLevelSystem(final ContinuousTPFESpace space)
			{
				final SparseMatrix s = new SparseMatrix(space.getShapeFunctions()
				                                             .size(),
				                                        space.getShapeFunctions()
				                                             .size());
				final DenseVector rhs = new DenseVector(space.getShapeFunctions()
				                                             .size());
				final TPCellIntegral<ContinuousTPShapeFunction> gg =
					new TPCellIntegral<>(ScalarFunction.constantFunction(1),
					                     TPCellIntegral.GRAD_GRAD);
				space.writeCellIntegralsToMatrix(List.of(gg), s);
				final TPRightHandSideIntegral<ContinuousTPShapeFunction> f =
					new TPRightHandSideIntegral<>(LaplaceReferenceSolution.scalarReferenceSolution(),
					                              TPRightHandSideIntegral.VALUE);
				space.writeCellIntegralsToRhs(List.of(f), rhs);
				space.writeBoundaryValuesTo(LaplaceReferenceSolution.scalarReferenceSolution(),
				                            s,
				                            rhs);
				return new Tuple2<>(s, rhs);
			}
			
			@Override
			public List<Smoother> createSmoothers()
			{
				final ArrayList<Smoother> ret = new ArrayList<>();
				for (int i = 1; i < spaces.size(); i++)
					ret.add(new RichardsonSmoother(0.1, 4));
				return ret;
			}
		};
		DenseVector solut = new DenseVector(mg.finest_rhs.mul(0));
		for (final IntCoordinates c : solut.getShape()
		                                   .range())
			solut.add(Math.random() * 0.1, c);
		PlotWindow.addPlot(new ScalarPlot2D(LaplaceReferenceSolution.scalarReferenceSolution(),
		                                    mg.spaces.get(refinements)
		                                             .generatePlotPoints(
			                                             30),
		                                    30, "reference"));
		ScalarFESpaceFunction<ContinuousTPShapeFunction> sol =
			new ScalarFESpaceFunction<>(
				mg.spaces.get(refinements)
				         .getShapeFunctions(), solut);
		PlotWindow.addPlot(new ScalarPlot2D(sol,
		                                    mg.spaces.get(refinements)
		                                             .generatePlotPoints(30),
		                                    30, "initial"));

//		solut = mg.systems.get(0)
//		                  .solve(mg.restrictionOperator.get(0)
//		                                               .mvMul(mg.restrictionOperator.get(1)
//		                                                                            .mvMul(mg.rhs)));
//		sol =
//			new ScalarFESpaceFunction<>(
//				mg.spaces.get(0)
//				         .getShapeFunctions(), solut);
//		PlotWindow.addPlot(new ScalarPlot2D(sol,
//		                                    mg.spaces.get(refinements)
//		                                             .generatePlotPoints(30),
//		                                    30, "coarse"));
//		solut = mg.systems.get(1)
//		                  .solve(mg.restrictionOperator.get(1)
//		                                               .mvMul(mg.rhs));
//		sol =
//			new ScalarFESpaceFunction<>(
//				mg.spaces.get(1)
//				         .getShapeFunctions(), solut);
//		PlotWindow.addPlot(new ScalarPlot2D(sol,
//		                                    mg.spaces.get(refinements)
//		                                             .generatePlotPoints(30),
//		                                    30, "coarse1"));
		//solut = new DenseVector(mg.rhs.mul(0));KC
//		final SparseMatrix A_0 = (SparseMatrix) mg.systems.get(0);
//		final SparseMatrix A_1 = (SparseMatrix) mg.systems.get(1);
//		final Matrix PTA_1P = mg.prolongationMatrices2.get(0)
//		                                              .tmMul(A_1)
//		                                              .mmMul(mg.prolongationMatrices2.get(0));
//		System.out.println(A_0.sub(PTA_1P)
//		                      .absMaxElement());
//		PlotWindow.addPlot(new MatrixPlot(A_0));
//		PlotWindow.addPlot(new MatrixPlot(PTA_1P));

//		for (int i = 0; i < 10; i++)
//		{
//			System.out.println(solut.getLength());
//			solut = new DenseVector(mg.vCycle(solut, mg.rhs));
//			sol =
//				new ScalarFESpaceFunction<>(
//					mg.spaces.get(refinements)
//					         .getShapeFunctions(), solut);
//			PlotWindow.addPlot(new ScalarPlot2D(sol,
//			                                    mg.spaces.get(refinements)
//			                                             .generatePlotPoints(30),
//			                                    30, "iteration" + i));
//		}
		final IterativeSolver it = new IterativeSolver();
		it.showProgress = true;
		solut = new DenseVector(it.solvePGMRES(mg.finest_system,
		                                       new MGPReconditioner(mg),
		                                       mg.finest_rhs,
		                                       1e-8));
		sol =
			new ScalarFESpaceFunction<>(
				mg.spaces.get(refinements)
				         .getShapeFunctions(), solut);
		PlotWindow.addPlot(new ScalarPlot2D(sol,
		                                    mg.spaces.get(refinements)
		                                             .generatePlotPoints(30),
		                                    30, "iteration"));
	}
}
