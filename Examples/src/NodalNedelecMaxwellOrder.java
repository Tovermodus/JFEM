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
import java.util.TreeMap;

public class NodalNedelecMaxwellOrder
{
	public static void main(String[] args)
	{
		
		PerformanceArguments.PerformanceArgumentBuilder builder =
			new PerformanceArguments.PerformanceArgumentBuilder();
		builder.build();CoordinateVector start = CoordinateVector.fromValues(0, 0,0);
		CoordinateVector end = CoordinateVector.fromValues(1, 1,1);
		int polynomialDegree = 2;
		TPVectorCellIntegral<NodalNedelecShapeFunction> valueValue =
			new TPVectorCellIntegral<>(TPVectorCellIntegral.ROT_ROT);
		MixedCellIntegral<TPCell,TPFace,TPEdge, ContinuousTPShapeFunction, NodalNedelecShapeFunction>
			divValue =
			new MixedTPCellIntegral<>(ScalarFunction.constantFunction(-1),
				MixedTPCellIntegral.DIV_VALUE);
		MixedCellIntegral<TPCell, TPFace, TPEdge,ContinuousTPShapeFunction, NodalNedelecShapeFunction> vv =
			MixedCellIntegral.fromVelocityIntegral(valueValue);
		List<CellIntegral<TPCell, MixedShapeFunction<TPCell, TPFace,TPEdge,ContinuousTPShapeFunction,
			NodalNedelecShapeFunction>>> cellIntegrals =
			new ArrayList<>();
		cellIntegrals.add(vv);
		cellIntegrals.add(divValue);
		List<FaceIntegral<TPFace, MixedShapeFunction<TPCell, TPFace,TPEdge,ContinuousTPShapeFunction,
			NodalNedelecShapeFunction>>> faceIntegrals = new ArrayList<>();
		MixedRightHandSideIntegral<TPCell, TPFace, TPEdge, ContinuousTPShapeFunction, NodalNedelecShapeFunction> rightHandSideIntegral =
			MixedRightHandSideIntegral.fromVelocityIntegral(
				new TPVectorRightHandSideIntegral<NodalNedelecShapeFunction>(MaxwellReferenceSolution.rightHandSide(),
					TPVectorRightHandSideIntegral.VALUE));
		List<RightHandSideIntegral<TPCell, MixedShapeFunction<TPCell, TPFace,TPEdge,
			ContinuousTPShapeFunction,
			NodalNedelecShapeFunction>>> rightHandSideIntegrals = new ArrayList<>();
		rightHandSideIntegrals.add(rightHandSideIntegral);
		
		List<BoundaryRightHandSideIntegral< TPFace, MixedShapeFunction<TPCell, TPFace,
			TPEdge,ContinuousTPShapeFunction,
			NodalNedelecShapeFunction>>> boundaryFaceIntegrals =
			new ArrayList<>();
		MixedNodalNedelecSpace grid = null;
		List<ScalarFunction> solutions = new ArrayList<>();
		List<VectorFunction> solutionsVec = new ArrayList<>();
		for(int i = 0; i < 3; i++)
		{
			grid = new MixedNodalNedelecSpace(start, end,
				Ints.asList(3*(int)Math.pow(2,i), 3*(int)Math.pow(2,i), 3*(int)Math.pow(2,i)),
				polynomialDegree);
			grid.assembleCells();
			grid.assembleFunctions(polynomialDegree);
			grid.initializeSystemMatrix();
			grid.initializeRhs();
			System.out.println("velocity Boundary");
			grid.setVelocityBoundaryValues(MaxwellReferenceSolution.vectorBoundaryValues());
			System.out.println("pressure Boundary");
			grid.setPressureBoundaryValues(MaxwellReferenceSolution.pressureBoundaryValues());
			System.out.println("Cell Integrals");
			grid.evaluateCellIntegrals(cellIntegrals, rightHandSideIntegrals);
			System.out.println("Face Integrals");
			grid.evaluateFaceIntegrals(faceIntegrals, boundaryFaceIntegrals);
			System.out.println("solve system: " + grid.getSystemMatrix().getRows() + "×" + grid.getSystemMatrix().getCols());
			IterativeSolver<SparseMatrix> it = new IterativeSolver<>();
			Vector solution1 = it.solveGMRES(grid.getSystemMatrix(), grid.getRhs(), 1e-7);
			MixedFESpaceFunction<TPCell,TPFace,TPEdge,ContinuousTPShapeFunction, NodalNedelecShapeFunction> solut =
				new MixedFESpaceFunction<>(
					grid.getShapeFunctions(), solution1);
			solutions.add(solut.getPressureFunction());
			solutionsVec.add(solut.getVelocityFunction());
		}
		solutions.add(MaxwellReferenceSolution.pressureReferenceSolution());
		solutionsVec.add(MaxwellReferenceSolution.velocityReferenceSolution());
		
		Map<String, Map<CoordinateVector, Double>> valList = new TreeMap<>();
		int k =0;
		for(ScalarFunction sol: solutions)
			valList.put("solution"+k++,sol.valuesInPoints(grid.generatePlotPoints(20)));
		for(VectorFunction sol: solutionsVec)
			valList.put("solutionVec"+k++,sol.componentValuesInPoints(grid.generatePlotPoints(20),0));
		System.out.println(ConvergenceOrderEstimator.estimateL2Vector(solutionsVec,
			grid.generatePlotPoints(20)));
		System.out.println(ConvergenceOrderEstimator.estimateL2Scalar(solutions, grid.generatePlotPoints(20)));
		new PlotFrame(valList,start, end);
	}
}
