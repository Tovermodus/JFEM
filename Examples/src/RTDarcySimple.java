import basic.*;
import com.google.common.primitives.Ints;
import linalg.CoordinateVector;
import linalg.IterativeSolver;
import linalg.MatrixPlot;
import linalg.Vector;
import mixed.*;
import tensorproduct.*;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

public class RTDarcySimple
{
	public static void main(String[] args)
	{
		
		PerformanceArguments.PerformanceArgumentBuilder builder =
			new PerformanceArguments.PerformanceArgumentBuilder();
		builder.build();CoordinateVector start = CoordinateVector.fromValues(0,0);
		CoordinateVector end = CoordinateVector.fromValues(1, 1);
		int polynomialDegree = 2;
		QkQkSpace grid = new QkQkSpace(start, end,
			Ints.asList(5,5));
		TPVectorCellIntegral<ContinuousTPVectorFunction> valueValue =
			new TPVectorCellIntegral<>(TPVectorCellIntegral.VALUE_VALUE);
		MixedCellIntegral<TPCell, ContinuousTPShapeFunction, ContinuousTPVectorFunction,QkQkFunction>
			divValue = new MixedTPCellIntegral<>(MixedTPCellIntegral.DIV_VALUE);
		MixedCellIntegral<TPCell, ContinuousTPShapeFunction, ContinuousTPVectorFunction,QkQkFunction> vv =
			MixedCellIntegral.fromVelocityIntegral(valueValue);
		List<CellIntegral<TPCell,  QkQkFunction>> cellIntegrals =
			new ArrayList<>();
		cellIntegrals.add(vv);
		cellIntegrals.add(divValue);
		List<FaceIntegral< TPFace, QkQkFunction>> faceIntegrals = new ArrayList<>();
		MixedRightHandSideIntegral<TPCell, ContinuousTPShapeFunction, ContinuousTPVectorFunction,QkQkFunction> rightHandSideIntegral =
			MixedRightHandSideIntegral.fromPressureIntegral(
				new TPRightHandSideIntegral<>(ScalarFunction.constantFunction(-1),
					TPRightHandSideIntegral.VALUE));
		List<RightHandSideIntegral<TPCell,  QkQkFunction>> rightHandSideIntegrals = new ArrayList<>();
		rightHandSideIntegrals.add(rightHandSideIntegral);
		List<BoundaryRightHandSideIntegral<TPFace, QkQkFunction>> boundaryFaceIntegrals =
			new ArrayList<>();
		
		MixedBoundaryRightHandSideIntegral<TPFace, ContinuousTPShapeFunction, ContinuousTPVectorFunction,QkQkFunction> dirichlet =
			MixedBoundaryRightHandSideIntegral.fromVelocityIntegral(
				new TPVectorBoundaryFaceIntegral<>(ScalarFunction.constantFunction(0),
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
		IterativeSolver i = new IterativeSolver();
		Vector solution1 = i.solveGMRES(grid.getSystemMatrix(), grid.getRhs(), 1e-10);
		MixedFESpaceFunction<QkQkFunction> solut =
			new MixedFESpaceFunction<>(
				grid.getShapeFunctions(), solution1);
		System.out.println(solution1.absMaxElement());
		
		PlotWindow p = new PlotWindow();
		p.addPlot(new MatrixPlot(grid.getSystemMatrix()));
		p.addPlot(new MixedPlot2D(solut, grid.generatePlotPoints(20),20));
		//for(QkQkFunction f:grid.getShapeFunctions().values())
		//	p.addPlot(new MixedPlot2D(f, grid.generatePlotPoints(50),50));
	}
}
