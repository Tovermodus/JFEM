import basic.*;
import com.google.common.primitives.Ints;
import linalg.CoordinateVector;
import linalg.IterativeSolver;
import linalg.Vector;
import mixed.*;
import tensorproduct.*;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

public class RTDarcy
{
	public static void main(final String[] args)
	{
		
		final PerformanceArguments.PerformanceArgumentBuilder builder =
			new PerformanceArguments.PerformanceArgumentBuilder();
		builder.build();
		final CoordinateVector start = CoordinateVector.fromValues(-1, -1, -1);
		final CoordinateVector end = CoordinateVector.fromValues(1, 1, 1);
		final int polynomialDegree = 2;
		final MixedRTSpace grid = new MixedRTSpace(start, end,
		                                           Ints.asList(3, 3, 3));
		final TPVectorCellIntegral<RTShapeFunction> valueValue =
			new TPVectorCellIntegral<>(TPVectorCellIntegral.VALUE_VALUE);
		final MixedCellIntegral<TPCell, ContinuousTPShapeFunction, RTShapeFunction, RTMixedFunction>
			divValue = new MixedTPCellIntegral<>(MixedTPCellIntegral.DIV_VALUE);
		final MixedCellIntegral<TPCell, ContinuousTPShapeFunction, RTShapeFunction, RTMixedFunction> vv =
			MixedCellIntegral.fromVelocityIntegral(valueValue);
		final List<CellIntegral<TPCell, RTMixedFunction>> cellIntegrals =
			new ArrayList<>();
		cellIntegrals.add(vv);
		cellIntegrals.add(divValue);
		final List<FaceIntegral<TPFace, RTMixedFunction>> faceIntegrals = new ArrayList<>();
		
		final MixedRightHandSideIntegral<TPCell, ContinuousTPShapeFunction, RTShapeFunction, RTMixedFunction>
			rightHandSideIntegral =
			MixedRightHandSideIntegral.fromPressureIntegral(
				new TPRightHandSideIntegral<ContinuousTPShapeFunction>(ScalarFunction.constantFunction(-1),
				                                                       TPRightHandSideIntegral.VALUE));
		
		final List<RightHandSideIntegral<TPCell, RTMixedFunction>> rightHandSideIntegrals = new ArrayList<>();
		rightHandSideIntegrals.add(rightHandSideIntegral);
		
		final List<BoundaryRightHandSideIntegral<TPFace, RTMixedFunction>> boundaryFaceIntegrals =
			new ArrayList<>();
		
		final MixedBoundaryRightHandSideIntegral<TPFace, ContinuousTPShapeFunction, RTShapeFunction, RTMixedFunction>
			dirichlet =
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
		System.out.println("solve system: " + grid.getSystemMatrix()
		                                          .getRows() + "×" + grid.getSystemMatrix()
		                                                                 .getCols());
		//grid.A.makeParallelReady(12);
		if (grid.getRhs()
		        .getLength() < 100)
		{
			System.out.println(grid.getSystemMatrix());
			System.out.println(grid.getRhs());
		}
		final IterativeSolver i = new IterativeSolver();
		final Vector solution1 = i.solveGMRES(grid.getSystemMatrix(), grid.getRhs(), 1e-10);
		//Vector solution1 = new DenseMatrix(grid.getSystemMatrix()).solve(grid.getRhs());
		System.out.println("solved");
		System.out.println(solution1);
		System.out.println(grid.getSystemMatrix()
		                       .mvMul(solution1));
		//grid.A.print_formatted();
		//grid.rhs.print_formatted();
		final MixedTPFESpaceFunction<RTMixedFunction> solut =
			new MixedTPFESpaceFunction<>(
				grid.getShapeFunctions(), solution1);
		final ArrayList<Map<CoordinateVector, Double>> valList = new ArrayList<>();
		/*valList.add(solut.pressureValuesInPoints(grid.generatePlotPoints(50)));
		valList.add(LaplaceReferenceSolution.scalarReferenceSolution().valuesInPoints(grid.generatePlotPoints(50)));
		valList.add(solut.velocityComponentsInPoints(grid.generatePlotPoints(50), 0));
		valList.add(LaplaceReferenceSolution.scalarReferenceSolution().getGradientFunction().componentValuesInPoints(grid.generatePlotPoints(50),0));
		valList.add(solut.velocityComponentsInPoints(grid.generatePlotPoints(50), 1));
		valList.add(LaplaceReferenceSolution.scalarReferenceSolution().getGradientFunction()
		.componentValuesInPoints(grid.generatePlotPoints(50),1));*/
		valList.add(solut.getPressureFunction()
		                 .valuesInPoints(grid.generatePlotPoints(20)));
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
