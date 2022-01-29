import basic.*;
import com.google.common.primitives.Ints;
import linalg.CoordinateVector;
import linalg.IterativeSolver;
import linalg.Vector;
import tensorproduct.*;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

public class VectorLaplace
{
	public static void main(final String[] args)
	{
		
		final PerformanceArguments.PerformanceArgumentBuilder builder =
			new PerformanceArguments.PerformanceArgumentBuilder();
		builder.build();
		final CoordinateVector start = CoordinateVector.fromValues(-1, -1);
		final CoordinateVector end = CoordinateVector.fromValues(1, 1);
		final int polynomialDegree = 2;
		final TPVectorFESpace grid = new TPVectorFESpace(start, end,
		                                                 Ints.asList(10, 10));
		final TPVectorCellIntegral<TPVectorFunction> gg =
			new TPVectorCellIntegral<>(TPVectorCellIntegral.GRAD_GRAD);
		final TPVectorFaceIntegral<TPVectorFunction> gj =
			new TPVectorFaceIntegral<>(ScalarFunction.constantFunction(1),
			                           TPVectorFaceIntegral.GRAD_AVERAGE_VALUE_NORMALAVERAGE);
		final TPVectorFaceIntegral<TPVectorFunction> jg =
			new TPVectorFaceIntegral<>(ScalarFunction.constantFunction(1),
			                           TPVectorFaceIntegral.VALUE_NORMALAVERAGE_GRAD_AVERAGE);
		final double penalty = 100000;
		final TPVectorFaceIntegral<TPVectorFunction> jj =
			new TPVectorFaceIntegral<>(ScalarFunction.constantFunction(penalty),
			                           TPVectorFaceIntegral.VALUE_NORMALAVERAGE_VALUE_NORMALAVERAGE);
		final ArrayList<CellIntegral<TPCell, TPVectorFunction>> cellIntegrals =
			new ArrayList<>();
		cellIntegrals.add(gg);
		final ArrayList<FaceIntegral<TPFace, TPVectorFunction>> faceIntegrals = new ArrayList<>();
		faceIntegrals.add(jj);
		faceIntegrals.add(gj);
		faceIntegrals.add(jg);
		final TPVectorRightHandSideIntegral<TPVectorFunction> rightHandSideIntegral =
			new TPVectorRightHandSideIntegral<>(StokesReferenceSolution.rightHandSide(),
			                                    TPVectorRightHandSideIntegral.VALUE);
		//new VectorFunction()
//			{
//				@Override
//				public int getDomainDimension()
//				{
//					return 2;
//				}
//
//				@Override
//				public CoordinateVector value(CoordinateVector pos)
//				{
//					return CoordinateVector.fromValues(4,4);
//				}
//			},
//				TPVectorRightHandSideIntegral.VALUE);
		final ArrayList<RightHandSideIntegral<TPCell, TPVectorFunction>> rightHandSideIntegrals = new ArrayList<>();
		rightHandSideIntegrals.add(rightHandSideIntegral);
		final ScalarFunction func = new ScalarFunction()
		{
			@Override
			public int getDomainDimension()
			{
				return 2;
			}
			
			@Override
			public Double value(final CoordinateVector pos)
			{
				if (Math.abs(pos.x()) == 1)
					return penalty * (2 - Math.max(1, 2 * Math.abs(pos.y())));
				if (Math.abs(pos.y()) == 1)
					return penalty * (2 - Math.max(1, 2 * Math.abs(pos.x())));
				return (double) 0;
			}
		};
		final TPVectorBoundaryFaceIntegral<TPVectorFunction> bound
			= new TPVectorBoundaryFaceIntegral<>(new VectorFunction()
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
				return CoordinateVector.fromValues(0, 0);//func.value(pos), func.value(pos));
			}
		}, TPVectorBoundaryFaceIntegral.VALUE);
		final ArrayList<BoundaryRightHandSideIntegral<TPFace, TPVectorFunction>> boundaryFaceIntegrals
			= new ArrayList<>();
		boundaryFaceIntegrals.add(bound);
		grid.assembleCells();
		grid.assembleFunctions(polynomialDegree);
		grid.initializeSystemMatrix();
		grid.initializeRhs();
		System.out.println("Cell Integrals");
		grid.evaluateCellIntegrals(cellIntegrals, rightHandSideIntegrals);
		System.out.println("Face Integrals");
		grid.evaluateFaceIntegrals(faceIntegrals, boundaryFaceIntegrals);
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
		final VectorFESpaceFunction<TPVectorFunction> solut =
			new VectorFESpaceFunction<>(
				grid.getShapeFunctionMap(), solution1);
		final Map<CoordinateVector, Double> vals =
			solut.componentValuesInPoints(grid.generatePlotPoints(50), 0);
		final Map<CoordinateVector, Double> vals2 =
			solut.componentValuesInPoints(grid.generatePlotPoints(50), 1);
		new PlotFrame(List.of(vals, vals2), start, end);
	}
}
