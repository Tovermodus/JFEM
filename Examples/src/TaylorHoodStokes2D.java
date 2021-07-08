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
import java.util.OptionalDouble;

public class TaylorHoodStokes2D
{
	public static void main(String[] args)
	{
		
		PerformanceArguments.PerformanceArgumentBuilder builder =
			new PerformanceArguments.PerformanceArgumentBuilder();
		builder.build();
		CoordinateVector start = CoordinateVector.fromValues(0, 0);
		CoordinateVector end = CoordinateVector.fromValues(1, 1);
		int polynomialDegree = 1;
		TaylorHoodSpace grid = new TaylorHoodSpace(start, end,
			Ints.asList(3,3), polynomialDegree);
		TPVectorCellIntegral<ContinuousTPVectorFunction> valueValue =
			new TPVectorCellIntegral<>(ScalarFunction.constantFunction(StokesReferenceSolution.reynolds),
				TPVectorCellIntegral.GRAD_GRAD);
		MixedCellIntegral<TPCell,TPFace,TPEdge,ContinuousTPShapeFunction, ContinuousTPVectorFunction>
			divValue =
			new MixedTPCellIntegral<>(ScalarFunction.constantFunction(-1),
				MixedTPCellIntegral.DIV_VALUE);
		MixedCellIntegral<TPCell, TPFace,TPEdge, ContinuousTPShapeFunction, ContinuousTPVectorFunction> vv =
			MixedCellIntegral.fromVelocityIntegral(valueValue);
		List<CellIntegral<TPCell,  MixedShapeFunction<TPCell, TPFace,TPEdge,ContinuousTPShapeFunction,
			ContinuousTPVectorFunction>>> cellIntegrals =
			new ArrayList<>();
		cellIntegrals.add(vv);
		cellIntegrals.add(divValue);
		List<FaceIntegral< TPFace, MixedShapeFunction<TPCell, TPFace,TPEdge,ContinuousTPShapeFunction,
					ContinuousTPVectorFunction>>> faceIntegrals = new ArrayList<>();
		//faceIntegrals.add(jj);
		//faceIntegrals.add(gj);
		//faceIntegrals.add(jg);
		MixedRightHandSideIntegral<TPCell, TPFace,TPEdge, ContinuousTPShapeFunction,
			ContinuousTPVectorFunction> rightHandSideIntegral =
			MixedRightHandSideIntegral.fromVelocityIntegral(
				new TPVectorRightHandSideIntegral<ContinuousTPVectorFunction>(ScalarFunction.constantFunction(1).makeIsotropicVectorFunction(),
					TPVectorRightHandSideIntegral.VALUE));
		List<RightHandSideIntegral<TPCell,  MixedShapeFunction<TPCell, TPFace,TPEdge,
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
		System.out.println("Cell Integrals");
		grid.evaluateCellIntegrals(cellIntegrals, rightHandSideIntegrals);
		System.out.println("Face Integrals");
		grid.evaluateFaceIntegrals(faceIntegrals, boundaryFaceIntegrals);
		grid.setVelocityBoundaryValues(ScalarFunction.constantFunction(0).makeIsotropicVectorFunction());
		//grid.setPressureBoundaryValues(ScalarFunction.constantFunction(0));
		//grid.A.makeParallelReady(12);
		
//		for(int i = 0; i < grid.getShapeFunctions().size(); i++)
//		{
//			grid.getSystemMatrix().set(0,0,i);
//		}
//		grid.getSystemMatrix().set(1,0,0);
//		grid.getRhs().set(0,0);
		if (grid.getRhs().getLength() < 100)
		{
			System.out.println(grid.getSystemMatrix());
			System.out.println(grid.getRhs());
			//throw new IllegalStateException();
		}
		System.out.println("solve system: " + grid.getSystemMatrix().getRows() + "×" + grid.getSystemMatrix().getCols());
		//DenseMatrix A = new DenseMatrix(grid.getSystemMatrix());
		IterativeSolver<SparseMatrix> i = new IterativeSolver<>();
		Vector solution1 = i.solveCG(grid.getSystemMatrix(), grid.getRhs(), 3.3e-6);
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
		PlotWindow p = new PlotWindow();
		p.addPlot(new MixedPlot2D(solut,grid.generatePlotPoints(20),20));
	}
}
