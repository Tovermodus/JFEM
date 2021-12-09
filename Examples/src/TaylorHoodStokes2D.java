import basic.PlotWindow;
import basic.ScalarFunction;
import basic.ScalarPlot2D;
import com.google.common.primitives.Ints;
import linalg.*;
import mixed.*;
import tensorproduct.ContinuousTPShapeFunction;
import tensorproduct.ContinuousTPVectorFunction;
import tensorproduct.TPVectorCellIntegral;
import tensorproduct.TPVectorRightHandSideIntegral;
import tensorproduct.geometry.TPCell;

import java.util.ArrayList;
import java.util.List;

public class TaylorHoodStokes2D
{
	public static void main(final String[] args)
	{
		final CoordinateVector start = CoordinateVector.fromValues(0, 0);
		final CoordinateVector end = CoordinateVector.fromValues(1, 1);
		final int polynomialDegree = 1;
		final TPVectorCellIntegral<ContinuousTPVectorFunction> gradGrad =
			new TPVectorCellIntegral<>(ScalarFunction.constantFunction(StokesReferenceSolution.reynolds),
			                           TPVectorCellIntegral.GRAD_GRAD);
		final MixedCellIntegral<TPCell, ContinuousTPShapeFunction, ContinuousTPVectorFunction, QkQkFunction>
			divValue =
			new MixedTPCellIntegral<>(ScalarFunction.constantFunction(-1),
			                          MixedTPCellIntegral.DIV_VALUE);
		final MixedCellIntegral<TPCell, ContinuousTPShapeFunction, ContinuousTPVectorFunction, QkQkFunction> vv
			=
			MixedCellIntegral.fromVelocityIntegral(gradGrad);
		final MixedRightHandSideIntegral<TPCell, ContinuousTPShapeFunction,
			ContinuousTPVectorFunction, QkQkFunction> rightHandSideIntegral =
			MixedRightHandSideIntegral.fromVelocityIntegral(
				new TPVectorRightHandSideIntegral<>(StokesReferenceSolution.rightHandSide(),
				                                    TPVectorRightHandSideIntegral.VALUE));
		
		final TaylorHoodSpace grid = new TaylorHoodSpace(start, end,
		                                                 Ints.asList(5, 5));
		grid.assembleCells();
		grid.assembleFunctions(polynomialDegree);
		grid.initializeSystemMatrix();
		grid.initializeRhs();
		System.out.println("Cell Integrals");
		grid.evaluateCellIntegrals(List.of(vv, divValue), List.of(rightHandSideIntegral));
		System.out.println("Face Integrals");
		grid.evaluateFaceIntegrals(new ArrayList<>(), new ArrayList<>());
		grid.setVelocityBoundaryValues(StokesReferenceSolution.vectorBoundaryValues());
		final SparseMatrix mat = grid.getSystemMatrix();
		final DenseVector rhs = grid.getRhs();
		IterativeSolver i = new IterativeSolver();
		i.showProgress = false;
		Vector solution1 = i.solveCG(grid.getSystemMatrix(), grid.getRhs(), 1e-10);
		
		MixedTPFESpaceFunction<QkQkFunction> solut =
			new MixedTPFESpaceFunction<>(
				grid.getShapeFunctions(), solution1);
		PlotWindow.addPlot(new MixedPlot2D(solut, grid.generatePlotPoints(20), 20));
		final TaylorHoodSpace grid2 = new TaylorHoodSpace(start, end,
		                                                  Ints.asList(5, 5));
		grid2.assembleCells();
		grid2.assembleFunctions(polynomialDegree);
		grid2.initializeSystemMatrix();
		grid2.initializeRhs();
		grid2.setVelocityBoundaryValues(StokesReferenceSolution.vectorBoundaryValues());
		System.out.println("Cell Integrals");
		grid2.evaluateCellIntegrals(List.of(vv, divValue), List.of(rightHandSideIntegral));
		System.out.println("Face Integrals");
		grid2.evaluateFaceIntegrals(new ArrayList<>(), new ArrayList<>());
		
		System.out.println(mat.sub(grid2.getSystemMatrix())
		                      .absMaxElement());
		System.out.println(rhs.sub(grid2.getRhs())
		                      .absMaxElement());
		
		PlotWindow.addPlot(new MatrixPlot(mat));
		PlotWindow.addPlot(new MatrixPlot(grid2.getSystemMatrix()));
		
		System.out.println("solve system: " + grid2.getSystemMatrix()
		                                           .getRows() + "×" + grid2.getSystemMatrix()
		                                                                   .getCols());
		i = new IterativeSolver();
		i.showProgress = false;
		solution1 = i.solveCG(grid2.getSystemMatrix(), grid2.getRhs(), 1e-10);
		System.out.println("solved");
		System.out.println(grid2.getSystemMatrix()
		                        .sub(grid2.getSystemMatrix()
		                                  .transpose())
		                        .absMaxElement());
		solut =
			new MixedTPFESpaceFunction<>(
				grid2.getShapeFunctions(), solution1);
		PlotWindow.addPlot(new MixedPlot2D(solut, grid2.generatePlotPoints(20), 20));
		PlotWindow.addPlot(new ScalarPlot2D(solut.getVelocityFunction()
		                                         .getDivergenceFunction(), grid2.generatePlotPoints(20), 20));
	}
}
