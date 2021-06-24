import basic.*;
import com.google.common.primitives.Ints;
import linalg.CoordinateVector;
import linalg.IterativeSolver;
import linalg.SparseMatrix;
import linalg.Vector;
import mixed.*;
import tensorproduct.*;

import java.util.*;

public class QkQkMaxwell
{
	public static void main(String[] args)
	{
		PerformanceArguments.createInstance(true,12,true);
		CoordinateVector start = CoordinateVector.fromValues(0, 0,0);
		CoordinateVector end = CoordinateVector.fromValues(1, 1,1);
		int polynomialDegree = 2;
		QkQkSpace grid = new QkQkSpace(start, end,
			Ints.asList(8,8,8), polynomialDegree);
		TPVectorCellIntegral<ContinuousTPVectorFunction> valueValue =
			new TPVectorCellIntegral<>(TPVectorCellIntegral.ROT_ROT);
		MixedCellIntegral<TPCell,TPFace,TPEdge,ContinuousTPShapeFunction, ContinuousTPVectorFunction>
			divValue =
			new MixedTPCellIntegral<>(ScalarFunction.constantFunction(-1),
				MixedTPCellIntegral.DIV_VALUE);
		MixedCellIntegral<TPCell, TPFace, TPEdge,ContinuousTPShapeFunction, ContinuousTPVectorFunction> vv =
			MixedCellIntegral.fromVelocityIntegral(valueValue);
		List<CellIntegral<TPCell,  MixedShapeFunction<TPCell, TPFace,TPEdge,ContinuousTPShapeFunction,
			ContinuousTPVectorFunction>>> cellIntegrals =
			new ArrayList<>();
		cellIntegrals.add(vv);
		cellIntegrals.add(divValue);
		List<FaceIntegral<TPFace, MixedShapeFunction<TPCell, TPFace,TPEdge,ContinuousTPShapeFunction,
			ContinuousTPVectorFunction>>> faceIntegrals = new ArrayList<>();
		MixedRightHandSideIntegral<TPCell, TPFace, TPEdge,ContinuousTPShapeFunction,
			ContinuousTPVectorFunction> rightHandSideIntegral =
			MixedRightHandSideIntegral.fromVelocityIntegral(
				new TPVectorRightHandSideIntegral<ContinuousTPVectorFunction>(MaxwellReferenceSolution.rightHandSide(),
					TPVectorRightHandSideIntegral.VALUE));
		List<RightHandSideIntegral<TPCell, MixedShapeFunction<TPCell, TPFace,TPEdge,
			ContinuousTPShapeFunction,
			ContinuousTPVectorFunction>>> rightHandSideIntegrals = new ArrayList<>();
		rightHandSideIntegrals.add(rightHandSideIntegral);
		
		List<BoundaryRightHandSideIntegral< TPFace, MixedShapeFunction<TPCell, TPFace,
			TPEdge,ContinuousTPShapeFunction,
			ContinuousTPVectorFunction>>> boundaryFaceIntegrals =
			new ArrayList<>();
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
		//grid.A.makeParallelReady(12);
		
		/*for(int i = 0; i < grid.getShapeFunctions().size(); i++)
		{
			grid.getSystemMatrix().set(0,0,i);
		}
		grid.getSystemMatrix().set(1,0,0);
		grid.getRhs().set(0,0);
		*/
		if (grid.getRhs().getLength() < 100)
		{
			System.out.println(grid.getSystemMatrix());
			System.out.println(grid.getRhs());
			//throw new IllegalStateException();
		}
		System.out.println("solve system: " + grid.getSystemMatrix().getRows() + "×" + grid.getSystemMatrix().getCols());
		//DenseMatrix A = new DenseMatrix(grid.getSystemMatrix());
		IterativeSolver<SparseMatrix> i = new IterativeSolver<>();
		Vector solution1 = i.solveGMRES(grid.getSystemMatrix(), grid.getRhs(), 1e-6);
		//Vector solution1 = grid.getSystemMatrix().solve(grid.getRhs());
		//Vector solution1 = new DenseMatrix(grid.getSystemMatrix()).solve(grid.getRhs());
		System.out.println("solved");
		System.out.println("sol"+solution1);
		System.out.println("rhs"+grid.getRhs());
		
		System.out.println("rhs2"+grid.getSystemMatrix().mvMul(solution1));
		//grid.A.print_formatted();
		//grid.rhs.print_formatted();
		MixedFESpaceFunction<TPCell,TPFace,TPEdge,ContinuousTPShapeFunction,ContinuousTPVectorFunction> solut =
			new MixedFESpaceFunction<>(
				grid.getShapeFunctions(), solution1);
		Map<String,Map<CoordinateVector, Double>> valList = new TreeMap<>();
//		valList.add(StokesReferenceSolution.pressureReferenceSolution().valuesInPoints(grid.generatePlotPoints(50)));
//		valList.add(StokesReferenceSolution.velocityReferenceSolution().componentValuesInPoints(grid.generatePlotPoints(50),0));
//		valList.add(StokesReferenceSolution.velocityReferenceSolution().componentValuesInPoints(grid.generatePlotPoints(50),1));
		valList.put("x-reference",
			MaxwellReferenceSolution.velocityReferenceSolution().componentValuesInPoints(grid.generatePlotPoints(20),0));
		valList.put("y-reference",
			MaxwellReferenceSolution.velocityReferenceSolution().componentValuesInPoints(grid.generatePlotPoints(20),1));
		valList.put("z-reference",
			MaxwellReferenceSolution.velocityReferenceSolution().componentValuesInPoints(grid.generatePlotPoints(20),2));
		valList.put("Pressure",solut.pressureValuesInPoints(grid.generatePlotPoints(20)));
		valList.put("x-velocity",solut.velocityComponentsInPoints(grid.generatePlotPoints(20), 0));
		valList.put("y-velocity",solut.velocityComponentsInPoints(grid.generatePlotPoints(20), 1));
		valList.put("z-velocity",solut.velocityComponentsInPoints(grid.generatePlotPoints(20), 2));
		
		/*for(MixedShapeFunction<TPCell, TPFace, ContinuousTPShapeFunction,ContinuousTPVectorFunction>
		shapeFunction:grid.getShapeFunctions().values())
		
		{
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
