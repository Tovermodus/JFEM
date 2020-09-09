import basic.*;
import com.google.common.primitives.Ints;
import linalg.CoordinateVector;
import linalg.IterativeSolver;
import linalg.SparseMatrix;
import linalg.Vector;
import tensorproduct.*;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

public class VectorLaplace
{
	public static void main(String[] args)
	{
		CoordinateVector start = CoordinateVector.fromValues(-1, -1);
		CoordinateVector end = CoordinateVector.fromValues(1, 1);
		int polynomialDegree = 3;
		TPVectorFESpace grid = new TPVectorFESpace(start, end,
			Ints.asList(10, 10), polynomialDegree);
		TPVectorCellIntegral<TPVectorFunction> gg =
			new TPVectorCellIntegral<>(TPVectorCellIntegral.GRAD_GRAD);
		TPVectorFaceIntegral<TPVectorFunction> gj =
			new TPVectorFaceIntegral<>(ScalarFunction.constantFunction(1),
				TPVectorFaceIntegral.GRAD_AVERAGE_VALUE_NORMALAVERAGE);
		TPVectorFaceIntegral<TPVectorFunction> jg =
			new TPVectorFaceIntegral<>(ScalarFunction.constantFunction(1),
				TPVectorFaceIntegral.VALUE_NORMALAVERAGE_GRAD_AVERAGE);
		TPVectorFaceIntegral<TPVectorFunction> jj =
			new TPVectorFaceIntegral<>(ScalarFunction.constantFunction(100000),
				TPVectorFaceIntegral.VALUE_NORMALAVERAGE_VALUE_NORMALAVERAGE);
		ArrayList<CellIntegral<TPCell, TPFace, TPVectorFunction>> cellIntegrals =
			new ArrayList<>();
		cellIntegrals.add(gg);
		ArrayList<FaceIntegral<TPCell, TPFace, TPVectorFunction>> faceIntegrals = new ArrayList<>();
		faceIntegrals.add(jj);
		faceIntegrals.add(gj);
		faceIntegrals.add(jg);
		TPVectorRightHandSideIntegral<TPVectorFunction> rightHandSideIntegral =
			new TPVectorRightHandSideIntegral<>(new VectorFunction()
			{
				@Override
				public int getDomainDimension()
				{
					return 2;
				}
				
				@Override
				public CoordinateVector value(CoordinateVector pos)
				{
					return CoordinateVector.fromValues(4,4);
				}
			},
				TPVectorRightHandSideIntegral.VALUE);
		ArrayList<RightHandSideIntegral<TPCell, TPFace, TPVectorFunction>> rightHandSideIntegrals = new ArrayList<>();
		rightHandSideIntegrals.add(rightHandSideIntegral);
		ArrayList<BoundaryRightHandSideIntegral<TPCell, TPFace, TPVectorFunction>> boundaryFaceIntegrals = new ArrayList<>();
		grid.assembleCells();
		grid.assembleFunctions(polynomialDegree);
		grid.initializeSystemMatrix();
		grid.initializeRhs();
		System.out.println("Cell Integrals");
		grid.evaluateCellIntegrals(cellIntegrals, rightHandSideIntegrals);
		System.out.println("Face Integrals");
		grid.evaluateFaceIntegrals(faceIntegrals, boundaryFaceIntegrals);
		System.out.println("solve system: " + grid.getSystemMatrix().getRows() + "×" + grid.getSystemMatrix().getCols());
		//grid.A.makeParallelReady(12);
		if (grid.getRhs().getLength() < 50)
		{
			System.out.println(grid.getSystemMatrix());
			System.out.println(grid.getRhs());
		}
		IterativeSolver<SparseMatrix> i = new IterativeSolver<>();
		System.out.println("start stopwatch");
		Vector solution1 = i.solveCG(grid.getSystemMatrix(), grid.getRhs(), 1e-3);
		//Vector solution = ((DenseMatrix)grid.getSystemMatrix()).solve(grid.getRhs());
		System.out.println("solved");
		//grid.A.print_formatted();
		//grid.rhs.print_formatted();
		VectorFESpaceFunction<TPVectorFunction> solut =
			new VectorFESpaceFunction<>(
				grid.getShapeFunctions(), solution1);
		Map<CoordinateVector, Double> vals =
			solut.componentValuesInPoints(grid.generatePlotPoints(50),0);
		new PlotFrame(List.of(vals),start,end);
	}
}
