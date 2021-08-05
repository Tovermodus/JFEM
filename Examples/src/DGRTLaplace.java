import basic.*;
import com.google.common.primitives.Ints;
import linalg.CoordinateVector;
import linalg.IterativeSolver;
import linalg.Vector;
import tensorproduct.*;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

import java.util.ArrayList;
import java.util.Map;

public class DGRTLaplace
{
	public static void main(String[] args)
	{
		PerformanceArguments.PerformanceArgumentBuilder builder =
			new PerformanceArguments.PerformanceArgumentBuilder();
		builder.build();
		CoordinateVector start = CoordinateVector.fromValues(-1, -1);
		CoordinateVector end = CoordinateVector.fromValues(1, 1);
		int polynomialDegree = 2;
		RTFESpace grid = new RTFESpace(start, end,
			Ints.asList(15, 15));
		TPVectorCellIntegral<RTShapeFunction> gg =
			new TPVectorCellIntegral<>(TPVectorCellIntegral.GRAD_GRAD);
		TPVectorFaceIntegral<RTShapeFunction> gj =
			new TPVectorFaceIntegral<>(ScalarFunction.constantFunction(1),
				TPVectorFaceIntegral.GRAD_AVERAGE_VALUE_NORMALAVERAGE);
		TPVectorFaceIntegral<RTShapeFunction> jg =
			new TPVectorFaceIntegral<>(ScalarFunction.constantFunction(1),
				TPVectorFaceIntegral.VALUE_NORMALAVERAGE_GRAD_AVERAGE);
		TPVectorFaceIntegral<RTShapeFunction> jj =
			new TPVectorFaceIntegral<>(ScalarFunction.constantFunction(100000),
				TPVectorFaceIntegral.VALUE_NORMALAVERAGE_VALUE_NORMALAVERAGE);
		ArrayList<CellIntegral<TPCell, RTShapeFunction>> cellIntegrals =
			new ArrayList<>();
		cellIntegrals.add(gg);
		ArrayList<FaceIntegral<TPFace, RTShapeFunction>> faceIntegrals = new ArrayList<>();
		//faceIntegrals.add(jj);
		//faceIntegrals.add(gj);
		//faceIntegrals.add(jg);
		TPVectorRightHandSideIntegral<RTShapeFunction> rightHandSideIntegral =
			new TPVectorRightHandSideIntegral<>(new VectorFunction()
			{
				@Override
				public int getRangeDimension()
				{
					return 2;
				}
				@Override
				public int getDomainDimension()
				{
					return 2;
				}
				
				@Override
				public CoordinateVector value(CoordinateVector pos)
				{
					return CoordinateVector.fromValues(4, 4);
				}
			},
				TPVectorRightHandSideIntegral.VALUE);
		ArrayList<RightHandSideIntegral<TPCell, RTShapeFunction>> rightHandSideIntegrals = new ArrayList<>();
		rightHandSideIntegrals.add(rightHandSideIntegral);
		ArrayList<BoundaryRightHandSideIntegral< TPFace, RTShapeFunction>> boundaryFaceIntegrals = new ArrayList<>();
		grid.assembleCells();
		grid.assembleFunctions(polynomialDegree);
		grid.initializeSystemMatrix();
		grid.initializeRhs();
		System.out.println("Cell Integrals");
		grid.evaluateCellIntegrals(cellIntegrals, rightHandSideIntegrals);
		System.out.println("Face Integrals");
		grid.evaluateFaceIntegrals(faceIntegrals, boundaryFaceIntegrals);
		grid.setBoundaryValues(new VectorFunction()
		{
			@Override
			public int getRangeDimension()
			{
				return 2;
			}
			@Override
			public int getDomainDimension()
			{
				return 2;
			}
			
			@Override
			public CoordinateVector value(CoordinateVector pos)
			{
				return CoordinateVector.fromValues(0,0);
			}
		});
		System.out.println("solve system: " + grid.getSystemMatrix().getRows() + "×" + grid.getSystemMatrix().getCols());
		//grid.A.makeParallelReady(12);
		if (grid.getRhs().getLength() < 50)
		{
			System.out.println(grid.getSystemMatrix());
			System.out.println(grid.getRhs());
		}
		IterativeSolver i = new IterativeSolver();
		System.out.println("start stopwatch");
		Vector solution1 = i.solveCG(grid.getSystemMatrix(), grid.getRhs(), 1e-3);
		//Vector solution = ((DenseMatrix)grid.getSystemMatrix()).solve(grid.getRhs());
		System.out.println("solved");
		//grid.A.print_formatted();
		//grid.rhs.print_formatted();
		VectorFESpaceFunction<RTShapeFunction> solut =
			new VectorFESpaceFunction<>(
				grid.getShapeFunctions(), solution1);
		ArrayList<Map<CoordinateVector, Double>> valList = new ArrayList<>();
		valList.add(solut.componentValuesInPoints(grid.generatePlotPoints(50), 0));
		valList.add(solut.componentValuesInPoints(grid.generatePlotPoints(50), 1));
		new PlotFrame(valList, start, end);
	}
}
