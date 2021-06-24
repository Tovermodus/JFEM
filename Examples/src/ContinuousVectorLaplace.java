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
		PerformanceArguments.createInstance(true,12,true);
		CoordinateVector start = CoordinateVector.fromValues(-1, -1);
		CoordinateVector end = CoordinateVector.fromValues(1, 1);
		int polynomialDegree = 2;
		ContinuousTPFEVectorSpace grid = new ContinuousTPFEVectorSpace(start, end,
			Ints.asList(20, 20), polynomialDegree);
		TPVectorCellIntegral<ContinuousTPVectorFunction> gg =
			new TPVectorCellIntegral<>(TPVectorCellIntegral.GRAD_GRAD);
		ArrayList<CellIntegral<TPCell, ContinuousTPVectorFunction>> cellIntegrals =
			new ArrayList<>();
		cellIntegrals.add(gg);
		TPVectorFaceIntegral<ContinuousTPVectorFunction> jj =
			new TPVectorFaceIntegral<>(ScalarFunction.constantFunction(1000),
			TPVectorFaceIntegral.BOUNDARY_VALUE);
		ArrayList<FaceIntegral<TPFace, ContinuousTPVectorFunction>> faceIntegrals = new ArrayList<>();
		//faceIntegrals.add(jj);
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
		ArrayList<RightHandSideIntegral<TPCell, ContinuousTPVectorFunction>> rightHandSideIntegrals = new ArrayList<>();
		rightHandSideIntegrals.add(rightHandSideIntegral);
		ArrayList<BoundaryRightHandSideIntegral<TPFace, ContinuousTPVectorFunction>> boundaryFaceIntegrals = new ArrayList<>();
		grid.assembleCells();
		grid.assembleFunctions(polynomialDegree);
		grid.initializeSystemMatrix();
		grid.initializeRhs();
		System.out.println("Cell Integrals");
		grid.evaluateCellIntegrals(cellIntegrals, rightHandSideIntegrals);
		System.out.println("Face Integrals");
		grid.evaluateFaceIntegrals(faceIntegrals, boundaryFaceIntegrals);
		System.out.println("Boundary Values");
		ScalarFunction func = new ScalarFunction()
		{
			@Override
			public int getDomainDimension()
			{
				return 2;
			}
			
			@Override
			public Double value(CoordinateVector pos)
			{
				if (Math.abs(pos.x()) == 1)
					return (2 - Math.max(1, 2 * Math.abs(pos.y())));
				if (Math.abs(pos.y()) == 1)
					return (2 - Math.max(1, 2 * Math.abs(pos.x())));
				return (double) 0;
			}
		};
		grid.setBoundaryValues(new VectorFunction()
		{
			@Override
			public int getDomainDimension()
			{
				return 2;
			}
			
			@Override
			public CoordinateVector value(CoordinateVector pos)
			{
				return CoordinateVector.fromValues(func.value(pos), func.value(pos));
			}
		});
		System.out.println("solve system: " + grid.getSystemMatrix().getRows() + "×" + grid.getSystemMatrix().getCols());
		//grid.A.makeParallelReady(12);
		if (grid.getRhs().getLength() < 50)
		{
			System.out.println(grid.getSystemMatrix());
			System.out.println(grid.getRhs());
		}
		IterativeSolver<SparseMatrix> i = new IterativeSolver<>();
		System.out.println("start stopwatch");
		Vector solution1 = i.solveGMRES(grid.getSystemMatrix(), grid.getRhs(), 1e-3);
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
