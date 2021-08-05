import basic.*;
import com.google.common.primitives.Ints;
import linalg.CoordinateVector;
import linalg.IterativeSolver;
import linalg.Vector;
import mixed.*;
import tensorproduct.*;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

import java.util.*;

public class QkQkMaxwell
{
	public static void main(String[] args)
	{
		
		PerformanceArguments.PerformanceArgumentBuilder builder =
			new PerformanceArguments.PerformanceArgumentBuilder();
		builder.build();CoordinateVector start = CoordinateVector.fromValues(0, 0,0);
		CoordinateVector end = CoordinateVector.fromValues(1, 1,1);
		int polynomialDegree = 2;
		QkQkSpace grid = new QkQkSpace(start, end,
			Ints.asList(8,8,8));
		TPVectorCellIntegral<ContinuousTPVectorFunction> valueValue =
			new TPVectorCellIntegral<>(TPVectorCellIntegral.ROT_ROT);
		MixedCellIntegral<TPCell, ContinuousTPShapeFunction, ContinuousTPVectorFunction,QkQkFunction>
			divValue =
			new MixedTPCellIntegral<>(ScalarFunction.constantFunction(-1),
				MixedTPCellIntegral.DIV_VALUE);
		MixedCellIntegral<TPCell, ContinuousTPShapeFunction, ContinuousTPVectorFunction,QkQkFunction> vv =
			MixedCellIntegral.fromVelocityIntegral(valueValue);
		List<CellIntegral<TPCell, QkQkFunction >> cellIntegrals =
			new ArrayList<>();
		cellIntegrals.add(vv);
		cellIntegrals.add(divValue);
		List<FaceIntegral<TPFace,QkQkFunction >> faceIntegrals = new ArrayList<>();
		MixedRightHandSideIntegral<TPCell, ContinuousTPShapeFunction,
			ContinuousTPVectorFunction,QkQkFunction> rightHandSideIntegral =
			MixedRightHandSideIntegral.fromVelocityIntegral(
				new TPVectorRightHandSideIntegral<ContinuousTPVectorFunction>(MaxwellReferenceSolution.rightHandSide(),
					TPVectorRightHandSideIntegral.VALUE));
		List<RightHandSideIntegral<TPCell, QkQkFunction>> rightHandSideIntegrals = new ArrayList<>();
		rightHandSideIntegrals.add(rightHandSideIntegral);
		
		List<BoundaryRightHandSideIntegral< TPFace, QkQkFunction>> boundaryFaceIntegrals =
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
		IterativeSolver i = new IterativeSolver();
		Vector solution1 = i.solveGMRES(grid.getSystemMatrix(), grid.getRhs(), 1e-6);
		//Vector solution1 = grid.getSystemMatrix().solve(grid.getRhs());
		//Vector solution1 = new DenseMatrix(grid.getSystemMatrix()).solve(grid.getRhs());
		System.out.println("solved");
		System.out.println("sol"+solution1);
		System.out.println("rhs"+grid.getRhs());
		
		System.out.println("rhs2"+grid.getSystemMatrix().mvMul(solution1));
		//grid.A.print_formatted();
		//grid.rhs.print_formatted();
		MixedFESpaceFunction<QkQkFunction> solut =
			new MixedFESpaceFunction<>(
				grid.getShapeFunctions(), solution1);
		
		PlotWindow p = new PlotWindow();
		p.addPlot(new MixedPlot3D(solut, grid.generatePlotPoints(20), 20));
		p.addPlot(new MixedPlot3D(MaxwellReferenceSolution.mixedReferenceSolution(), grid.generatePlotPoints(20), 20));
	}
}
