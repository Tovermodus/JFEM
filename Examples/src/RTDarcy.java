import basic.*;
import com.google.common.primitives.Ints;
import linalg.CoordinateVector;
import linalg.IterativeSolver;
import linalg.SparseMatrix;
import linalg.Vector;
import mixed.*;
import tensorproduct.*;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

public class RTDarcy
{
	public static void main(String[] args)
	{
		PerformanceArguments.createInstance(true,12,true);
		CoordinateVector start = CoordinateVector.fromValues(-1, -1,-1);
		CoordinateVector end = CoordinateVector.fromValues(1, 1,1);
		int polynomialDegree = 2;
		MixedRTSpace grid = new MixedRTSpace(start, end,
			Ints.asList(3,3,3), polynomialDegree);
		TPVectorCellIntegral<RTShapeFunction> valueValue =
			new TPVectorCellIntegral<>(TPVectorCellIntegral.VALUE_VALUE);
		MixedCellIntegral<TPCell,TPFace,TPEdge,ContinuousTPShapeFunction, RTShapeFunction>
			divValue = new MixedTPCellIntegral<>(MixedTPCellIntegral.DIV_VALUE);
		MixedCellIntegral<TPCell, TPFace,TPEdge, ContinuousTPShapeFunction, RTShapeFunction> vv =
			MixedCellIntegral.fromVelocityIntegral(valueValue);
		List<CellIntegral<TPCell,  MixedShapeFunction<TPCell, TPFace,TPEdge,ContinuousTPShapeFunction,
			RTShapeFunction>>> cellIntegrals =
			new ArrayList<>();
		cellIntegrals.add(vv);
		cellIntegrals.add(divValue);
		List<FaceIntegral< TPFace, MixedShapeFunction<TPCell, TPFace,TPEdge,ContinuousTPShapeFunction,
			RTShapeFunction>>> faceIntegrals = new ArrayList<>();
		
		MixedRightHandSideIntegral<TPCell, TPFace,TPEdge, ContinuousTPShapeFunction, RTShapeFunction> rightHandSideIntegral =
			MixedRightHandSideIntegral.fromPressureIntegral(
				new TPRightHandSideIntegral<ContinuousTPShapeFunction>(ScalarFunction.constantFunction(-1),
					TPRightHandSideIntegral.VALUE, false));
		
		
		List<RightHandSideIntegral<TPCell,  MixedShapeFunction<TPCell, TPFace,TPEdge,
			ContinuousTPShapeFunction,
			RTShapeFunction>>> rightHandSideIntegrals = new ArrayList<>();
		rightHandSideIntegrals.add(rightHandSideIntegral);
		
		
		List<BoundaryRightHandSideIntegral< TPFace, MixedShapeFunction<TPCell, TPFace,
			TPEdge,ContinuousTPShapeFunction,
			RTShapeFunction>>> boundaryFaceIntegrals =
			new ArrayList<>();
		
		MixedBoundaryRightHandSideIntegral<TPCell, TPFace, TPEdge,ContinuousTPShapeFunction, RTShapeFunction> dirichlet =
			MixedBoundaryRightHandSideIntegral.fromVelocityIntegral(
				new TPVectorBoundaryFaceIntegral<RTShapeFunction>(ScalarFunction.constantFunction(0),
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
		System.out.println("solve system: " + grid.getSystemMatrix().getRows() + "×" + grid.getSystemMatrix().getCols());
		//grid.A.makeParallelReady(12);
		if (grid.getRhs().getLength() < 100)
		{
			System.out.println(grid.getSystemMatrix());
			System.out.println(grid.getRhs());
		}
		IterativeSolver<SparseMatrix> i = new IterativeSolver<>();
		Vector solution1 = i.solveGMRES(grid.getSystemMatrix(), grid.getRhs(), 1e-10);
		//Vector solution1 = new DenseMatrix(grid.getSystemMatrix()).solve(grid.getRhs());
		System.out.println("solved");
		System.out.println(solution1);
		System.out.println(grid.getSystemMatrix().mvMul(solution1));
		//grid.A.print_formatted();
		//grid.rhs.print_formatted();
		MixedFESpaceFunction<TPCell,TPFace,TPEdge,ContinuousTPShapeFunction,RTShapeFunction> solut =
			new MixedFESpaceFunction<>(
				grid.getShapeFunctions(), solution1);
		ArrayList<Map<CoordinateVector, Double>> valList = new ArrayList<>();
		/*valList.add(solut.pressureValuesInPoints(grid.generatePlotPoints(50)));
		valList.add(LaplaceReferenceSolution.scalarReferenceSolution().valuesInPoints(grid.generatePlotPoints(50)));
		valList.add(solut.velocityComponentsInPoints(grid.generatePlotPoints(50), 0));
		valList.add(LaplaceReferenceSolution.scalarReferenceSolution().getGradientFunction().componentValuesInPoints(grid.generatePlotPoints(50),0));
		valList.add(solut.velocityComponentsInPoints(grid.generatePlotPoints(50), 1));
		valList.add(LaplaceReferenceSolution.scalarReferenceSolution().getGradientFunction()
		.componentValuesInPoints(grid.generatePlotPoints(50),1));*/
		valList.add(solut.getPressureFunction().valuesInPoints(grid.generatePlotPoints(20)));
		//valList.add(solut.getVelocityFunction().getDivergenceFunction().valuesInPoints(grid
		// .generatePlotPoints(20)));
		//valList.add(LaplaceReferenceSolution.scalarReferenceSolution().getGradientFunction()
		// .getDivergenceFunction().valuesInPoints(grid.generatePlotPoints(20)));
		/*for(int k = 0; k < grid.getShapeFunctions().size(); k++)
		{
			MixedShapeFunction<TPCell,TPFace,ContinuousTPShapeFunction,RTShapeFunction> shapeFunction =
				grid.getShapeFunctions().get(k);
			if(shapeFunction.isPressure())
			valList.add(shapeFunction.pressureValuesInPoints(grid.generatePlotPoints(50)));
			if(shapeFunction.isVelocity())
			{
				valList.add(shapeFunction.velocityComponentsInPoints(grid.generatePlotPoints(50),
					shapeFunction.getVelocityShapeFunction().getComponent()));
			}
		}*/
		new PlotFrame(valList, start, end);
	}
}
