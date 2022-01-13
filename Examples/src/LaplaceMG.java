import basic.PlotWindow;
import basic.ScalarFESpaceFunction;
import basic.ScalarPlot2D;
import io.vavr.Tuple2;
import linalg.*;
import multigrid.RichardsonSmoother;
import multigrid.Smoother;
import tensorproduct.*;

import java.util.ArrayList;
import java.util.List;

public class LaplaceMG
{
	public static void main(final String[] args)
	{
		final int refinements = 6;
		final TPCellIntegral<ContinuousTPShapeFunction> gg =
			new TPCellIntegral<>(TPCellIntegral.GRAD_GRAD);
		final TPRightHandSideIntegral<ContinuousTPShapeFunction> rightHandSideIntegral =
			new TPRightHandSideIntegral<>(LaplaceReferenceSolution.scalarRightHandSide(),
			                              TPRightHandSideIntegral.VALUE);
		final TPMultiGridSpace<ContinuousTPFESpace, ContinuousTPShapeFunction, Double, CoordinateVector, CoordinateMatrix>
			mg = new TPMultiGridSpace<>(CoordinateVector.fromValues(-1, -1),
			                            CoordinateVector.fromValues(1, 1),
			                            new IntCoordinates(3, 3),
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
			public Tuple2<VectorMultiplyable, DenseVector> createSystem(final ContinuousTPFESpace space)
			{
				final SparseMatrix s = new SparseMatrix(space.getShapeFunctions()
				                                             .size(),
				                                        space.getShapeFunctions()
				                                             .size());
				final DenseVector rhs = new DenseVector(space.getShapeFunctions()
				                                             .size());
				space.writeCellIntegralsToMatrix(List.of(gg), s);
				space.writeCellIntegralsToRhs(List.of(rightHandSideIntegral), rhs);
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
		
		final IterativeSolver it = new IterativeSolver();
		it.showProgress = true;
		solut = new DenseVector(it.solvePGMRES(mg.finest_system,
		                                       new TPMGPreconditioner(mg),
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
		System.out.println(solut.getLength());
	}
}
