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
	public static void main(final String[] args)
	{
		final PerformanceArguments.PerformanceArgumentBuilder builder =
			new PerformanceArguments.PerformanceArgumentBuilder();
		builder.build();
		final CoordinateVector start = CoordinateVector.fromValues(-1, -1);
		final CoordinateVector end = CoordinateVector.fromValues(1, 1);
		final int polynomialDegree = 2;
		final RTFESpace grid = new RTFESpace(start, end,
		                                     Ints.asList(15, 15));
		final TPVectorCellIntegral<RTShapeFunction> gg =
			new TPVectorCellIntegral<>(TPVectorCellIntegral.GRAD_GRAD);
		final TPVectorFaceIntegral<RTShapeFunction> gj =
			new TPVectorFaceIntegral<>(ScalarFunction.constantFunction(1),
			                           TPVectorFaceIntegral.GRAD_AVERAGE_VALUE_NORMALAVERAGE);
		final TPVectorFaceIntegral<RTShapeFunction> jg =
			new TPVectorFaceIntegral<>(ScalarFunction.constantFunction(1),
			                           TPVectorFaceIntegral.VALUE_NORMALAVERAGE_GRAD_AVERAGE);
		final TPVectorFaceIntegral<RTShapeFunction> jj =
			new TPVectorFaceIntegral<>(ScalarFunction.constantFunction(100000),
			                           TPVectorFaceIntegral.VALUE_NORMALAVERAGE_VALUE_NORMALAVERAGE);
		final ArrayList<CellIntegral<TPCell, RTShapeFunction>> cellIntegrals =
			new ArrayList<>();
		cellIntegrals.add(gg);
		final ArrayList<FaceIntegral<TPFace, RTShapeFunction>> faceIntegrals = new ArrayList<>();
		//faceIntegrals.add(jj);
		//faceIntegrals.add(gj);
		//faceIntegrals.add(jg);
		final TPVectorRightHandSideIntegral<RTShapeFunction> rightHandSideIntegral =
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
				public CoordinateVector value(final CoordinateVector pos)
				{
					return CoordinateVector.fromValues(4, 4);
				}
			},
			                                    TPVectorRightHandSideIntegral.VALUE);
		final ArrayList<RightHandSideIntegral<TPCell, RTShapeFunction>> rightHandSideIntegrals = new ArrayList<>();
		rightHandSideIntegrals.add(rightHandSideIntegral);
		final ArrayList<BoundaryRightHandSideIntegral<TPFace, RTShapeFunction>> boundaryFaceIntegrals
			= new ArrayList<>();
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
			public CoordinateVector value(final CoordinateVector pos)
			{
				return CoordinateVector.fromValues(0, 0);
			}
		});
		System.out.println("solve system: " + grid.getSystemMatrix()
		                                          .getRows() + "×" + grid.getSystemMatrix()
		                                                                 .getCols());
		//grid.A.makeParallelReady(12);
		if (grid.getRhs()
		        .getLength() < 50)
		{
			System.out.println(grid.getSystemMatrix());
			System.out.println(grid.getRhs());
		}
		final IterativeSolver i = new IterativeSolver();
		System.out.println("start stopwatch");
		final Vector solution1 = i.solveCG(grid.getSystemMatrix(), grid.getRhs(), 1e-3);
		//Vector solution = ((DenseMatrix)grid.getSystemMatrix()).solve(grid.getRhs());
		System.out.println("solved");
		//grid.A.print_formatted();
		//grid.rhs.print_formatted();
		final VectorFESpaceFunction<RTShapeFunction> solut =
			new VectorFESpaceFunction<>(
				grid.getShapeFunctionMap(), solution1);
		final ArrayList<Map<CoordinateVector, Double>> valList = new ArrayList<>();
		valList.add(solut.componentValuesInPoints(grid.generatePlotPoints(50), 0));
		valList.add(solut.componentValuesInPoints(grid.generatePlotPoints(50), 1));
		new PlotFrame(valList, start, end);
	}
}
