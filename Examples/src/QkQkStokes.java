import basic.*;
import com.google.common.primitives.Ints;
import linalg.CoordinateVector;
import linalg.IterativeSolver;
import linalg.Vector;
import mixed.*;
import tensorproduct.ContinuousTPShapeFunction;
import tensorproduct.ContinuousTPVectorFunction;
import tensorproduct.TPVectorCellIntegral;
import tensorproduct.TPVectorRightHandSideIntegral;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

public class QkQkStokes
{
	public static void main(final String[] args)
	{
		
		final PerformanceArguments.PerformanceArgumentBuilder builder =
			new PerformanceArguments.PerformanceArgumentBuilder();
		builder.build();
		final CoordinateVector start = CoordinateVector.fromValues(0, 0);
		final CoordinateVector end = CoordinateVector.fromValues(1, 1);
		final int polynomialDegree = 3;
		final QkQkSpace grid = new QkQkSpace(start, end,
		                                     Ints.asList(20, 20));
		final TPVectorCellIntegral<ContinuousTPVectorFunction> valueValue =
			new TPVectorCellIntegral<>(ScalarFunction.constantFunction(StokesReferenceSolution.reynolds),
			                           TPVectorCellIntegral.GRAD_GRAD);
		final MixedCellIntegral<TPCell, ContinuousTPShapeFunction, ContinuousTPVectorFunction, QkQkFunction>
			divValue =
			new MixedTPCellIntegral<>(ScalarFunction.constantFunction(-1),
			                          MixedTPCellIntegral.DIV_VALUE);
		final MixedCellIntegral<TPCell, ContinuousTPShapeFunction, ContinuousTPVectorFunction, QkQkFunction> vv
			=
			MixedCellIntegral.fromVelocityIntegral(valueValue);
		final List<CellIntegral<TPCell, QkQkFunction>> cellIntegrals =
			new ArrayList<>();
		cellIntegrals.add(vv);
		cellIntegrals.add(divValue);
		final List<FaceIntegral<TPFace, QkQkFunction>> faceIntegrals = new ArrayList<>();
		//faceIntegrals.add(jj);
		//faceIntegrals.add(gj);
		//faceIntegrals.add(jg);
		final MixedRightHandSideIntegral<TPCell, ContinuousTPShapeFunction,
			ContinuousTPVectorFunction, QkQkFunction> rightHandSideIntegral =
			MixedRightHandSideIntegral.fromVelocityIntegral(
				new TPVectorRightHandSideIntegral<ContinuousTPVectorFunction>(StokesReferenceSolution.rightHandSide(),
				                                                              TPVectorRightHandSideIntegral.VALUE));
		final List<RightHandSideIntegral<TPCell, QkQkFunction>> rightHandSideIntegrals = new ArrayList<>();
		rightHandSideIntegrals.add(rightHandSideIntegral);
		
		final List<BoundaryRightHandSideIntegral<TPFace, QkQkFunction>> boundaryFaceIntegrals =
			new ArrayList<>();
		grid.assembleCells();
		grid.assembleFunctions(polynomialDegree);
		grid.initializeSystemMatrix();
		grid.initializeRhs();
		System.out.println("Cell Integrals");
		grid.evaluateCellIntegrals(cellIntegrals, rightHandSideIntegrals);
		System.out.println("Face Integrals");
		grid.evaluateFaceIntegrals(faceIntegrals, boundaryFaceIntegrals);
		grid.setVelocityBoundaryValues(StokesReferenceSolution.vectorBoundaryValues());
		//grid.setPressureBoundaryValues(ScalarFunction.constantFunction(0));
		//grid.A.makeParallelReady(12);
		
		/*for(int i = 0; i < grid.getShapeFunctions().size(); i++)
		{
			grid.getSystemMatrix().set(0,0,i);
		}
		grid.getSystemMatrix().set(1,0,0);
		grid.getRhs().set(0,0);
		*/
		if (grid.getRhs()
		        .getLength() < 100)
		{
			System.out.println(grid.getSystemMatrix());
			System.out.println(grid.getRhs());
			//throw new IllegalStateException();
		}
		System.out.println("solve system: " + grid.getSystemMatrix()
		                                          .getRows() + "×" + grid.getSystemMatrix()
		                                                                 .getCols());
		//DenseMatrix A = new DenseMatrix(grid.getSystemMatrix());
		final IterativeSolver i = new IterativeSolver();
		final Vector solution1 = i.solveGMRES(grid.getSystemMatrix(), grid.getRhs(), 1e-7);
		//Vector solution1 = grid.getSystemMatrix().solve(grid.getRhs());
		//Vector solution1 = new DenseMatrix(grid.getSystemMatrix()).solve(grid.getRhs());
		System.out.println("solved");
		System.out.println("sol" + solution1);
		System.out.println("rhs" + grid.getRhs());
		
		System.out.println("rhs2" + grid.getSystemMatrix()
		                                .mvMul(solution1));
		//grid.A.print_formatted();
		//grid.rhs.print_formatted();
		final MixedTPFESpaceFunction<QkQkFunction> solut =
			new MixedTPFESpaceFunction<>(
				grid.getShapeFunctionMap(), solution1);
		final ArrayList<Map<CoordinateVector, Double>> valList = new ArrayList<>();
//		valList.add(StokesReferenceSolution.pressureReferenceSolution().valuesInPoints(grid.generatePlotPoints(50)));
//		valList.add(StokesReferenceSolution.velocityReferenceSolution().componentValuesInPoints(grid.generatePlotPoints(50),0));
//		valList.add(StokesReferenceSolution.velocityReferenceSolution().componentValuesInPoints(grid.generatePlotPoints(50),1));
		valList.add(solut.pressureValuesInPoints(grid.generatePlotPoints(50)));
		valList.add(solut.velocityComponentsInPoints(grid.generatePlotPoints(50), 0));
		valList.add(solut.velocityComponentsInPoints(grid.generatePlotPoints(50), 1));
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
