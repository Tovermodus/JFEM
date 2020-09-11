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

public class ContinuousVectorLaplace
{
	public static void main(String[] args)
	{
		CoordinateVector start = CoordinateVector.fromValues(-1, -1);
		CoordinateVector end = CoordinateVector.fromValues(1, 1);
		int polynomialDegree = 2;
		ContinuousTPFEVectorSpace grid = new ContinuousTPFEVectorSpace(start, end,
			Ints.asList(20, 20), polynomialDegree);
		TPVectorCellIntegral<ContinuousTPVectorFunction> gg =
			new TPVectorCellIntegral<>(TPVectorCellIntegral.GRAD_GRAD);
		ArrayList<CellIntegral<TPCell, TPFace, ContinuousTPVectorFunction>> cellIntegrals =
			new ArrayList<>();
		cellIntegrals.add(gg);
		TPVectorFaceIntegral<ContinuousTPVectorFunction> jj =
			new TPVectorFaceIntegral<>(ScalarFunction.constantFunction(1000),
			TPVectorFaceIntegral.BOUNDARY_VALUE);
		ArrayList<FaceIntegral<TPCell, TPFace, ContinuousTPVectorFunction>> faceIntegrals = new ArrayList<>();
		faceIntegrals.add(jj);
		TPVectorRightHandSideIntegral<ContinuousTPVectorFunction> rightHandSideIntegral =
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
		ArrayList<RightHandSideIntegral<TPCell, TPFace, ContinuousTPVectorFunction>> rightHandSideIntegrals = new ArrayList<>();
		rightHandSideIntegrals.add(rightHandSideIntegral);
		ArrayList<BoundaryRightHandSideIntegral<TPCell, TPFace, ContinuousTPVectorFunction>> boundaryFaceIntegrals = new ArrayList<>();
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
		VectorFESpaceFunction<ContinuousTPVectorFunction> solut =
			new VectorFESpaceFunction<>(
				grid.getShapeFunctions(), solution1);
		Map<CoordinateVector, Double> vals =
			solut.componentValuesInPoints(grid.generatePlotPoints(50),0);
		new PlotFrame(List.of(vals),start,end);
	}
}
