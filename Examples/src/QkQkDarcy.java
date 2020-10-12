import basic.*;
import com.google.common.primitives.Ints;
import linalg.*;
import mixed.*;
import tensorproduct.*;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

public class QkQkDarcy
{
	public static void main(String[] args)
	{
		CoordinateVector start = CoordinateVector.fromValues(-1, -1);
		CoordinateVector end = CoordinateVector.fromValues(1, 1);
		int polynomialDegree = 1;
		QkQkSpace grid = new QkQkSpace(start, end,
			Ints.asList(10, 10), polynomialDegree);TPVectorCellIntegral<ContinuousTPVectorFunction> valueValue =
		new TPVectorCellIntegral<>(TPVectorCellIntegral.VALUE_VALUE);
		MixedCellIntegral<TPCell, TPFace,TPEdge,ContinuousTPShapeFunction, ContinuousTPVectorFunction>
			divValue = new MixedTPCellIntegral<>(MixedTPCellIntegral.DIV_VALUE);
		MixedCellIntegral<TPCell, TPFace,TPEdge, ContinuousTPShapeFunction, ContinuousTPVectorFunction> vv =
			MixedCellIntegral.fromVelocityIntegral(valueValue);
		List<CellIntegral<TPCell, MixedShapeFunction<TPCell, TPFace,TPEdge,ContinuousTPShapeFunction,
			ContinuousTPVectorFunction>>> cellIntegrals =
			new ArrayList<>();
		cellIntegrals.add(vv);
		cellIntegrals.add(divValue);
		List<FaceIntegral< TPFace, MixedShapeFunction<TPCell, TPFace,TPEdge,ContinuousTPShapeFunction,
			ContinuousTPVectorFunction>>> faceIntegrals = new ArrayList<>();
		MixedRightHandSideIntegral<TPCell, TPFace,TPEdge, ContinuousTPShapeFunction,
			ContinuousTPVectorFunction> rightHandSideIntegral =
			MixedRightHandSideIntegral.fromPressureIntegral(
				new TPRightHandSideIntegral<ContinuousTPShapeFunction>(
					LaplaceReferenceSolution.scalarRightHandSide(),TPRightHandSideIntegral.VALUE, false));
		List<RightHandSideIntegral<TPCell, MixedShapeFunction<TPCell, TPFace,TPEdge,ContinuousTPShapeFunction,
			ContinuousTPVectorFunction>>> rightHandSideIntegrals = new ArrayList<>();
		rightHandSideIntegrals.add(rightHandSideIntegral);
		List<BoundaryRightHandSideIntegral< TPFace, MixedShapeFunction<TPCell, TPFace,
			TPEdge,ContinuousTPShapeFunction,
			ContinuousTPVectorFunction>>> boundaryFaceIntegrals = new ArrayList<>();
		MixedBoundaryRightHandSideIntegral<TPCell, TPFace, TPEdge,ContinuousTPShapeFunction,
			ContinuousTPVectorFunction> dirichlet =
			MixedBoundaryRightHandSideIntegral.fromVelocityIntegral(
				new TPVectorBoundaryFaceIntegral<ContinuousTPVectorFunction>(LaplaceReferenceSolution.scalarBoundaryValues(),
					TPVectorBoundaryFaceIntegral.NORMAL_VALUE));
		boundaryFaceIntegrals.add(dirichlet);
		grid.assembleCells();
		grid.assembleFunctions(polynomialDegree);
		grid.initializeSystemMatrix();
		grid.initializeRhs();
		System.out.println("Cell Integrals");
		grid.evaluateCellIntegrals(cellIntegrals, rightHandSideIntegrals);
		System.out.println("Face Integrals");
		grid.evaluateFaceIntegrals(faceIntegrals, boundaryFaceIntegrals);
		grid.setPressureBoundaryValues(ScalarFunction.constantFunction(0));
		System.out.println("solve system: " + grid.getSystemMatrix().getRows() + "×" + grid.getSystemMatrix().getCols());
		//grid.A.makeParallelReady(12);
		if (grid.getRhs().getLength() < 100)
		{
			System.out.println(grid.getSystemMatrix());
			System.out.println(grid.getRhs());
		}
		IterativeSolver<SparseMatrix> i = new IterativeSolver<>();
		Vector solution1 = i.solveGMRES(grid.getSystemMatrix(), grid.getRhs(), 1e-3);
		//Vector solution1 = new DenseMatrix(grid.getSystemMatrix()).solve(grid.getRhs());
		System.out.println("solved");
		System.out.println(solution1);
		System.out.println(grid.getSystemMatrix().mvMul(solution1));
		//grid.A.print_formatted();
		//grid.rhs.print_formatted();
		MixedFESpaceFunction<TPCell,TPFace,TPEdge,ContinuousTPShapeFunction,ContinuousTPVectorFunction> solut =
			new MixedFESpaceFunction<>(
				grid.getShapeFunctions(), solution1);
		ArrayList<Map<CoordinateVector, Double>> valList = new ArrayList<>();
		valList.add(solut.pressureValuesInPoints(grid.generatePlotPoints(50)));
		valList.add(solut.velocityComponentsInPoints(grid.generatePlotPoints(50), 0));
		valList.add(solut.velocityComponentsInPoints(grid.generatePlotPoints(50), 1));
		for(MixedShapeFunction<TPCell, TPFace,TPEdge, ContinuousTPShapeFunction,ContinuousTPVectorFunction> shapeFunction:grid.getShapeFunctions().values())
		{
			if(shapeFunction.isPressure())
			valList.add(shapeFunction.pressureValuesInPoints(grid.generatePlotPoints(50)));
			if(shapeFunction.isVelocity())
			{
				valList.add(shapeFunction.velocityComponentsInPoints(grid.generatePlotPoints(50),
					shapeFunction.getVelocityShapeFunction().getComponent()));
			}
		}
		new PlotFrame(valList, start, end);
	}
}
