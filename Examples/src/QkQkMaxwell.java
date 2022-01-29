import basic.*;
import com.google.common.primitives.Ints;
import linalg.BlockSparseMatrix;
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

public class QkQkMaxwell
{
	public static void main(final String[] args)
	{
		
		final PerformanceArguments.PerformanceArgumentBuilder builder =
			new PerformanceArguments.PerformanceArgumentBuilder();
		builder.build();
		final CoordinateVector start = CoordinateVector.fromValues(0, 0, 0);
		final CoordinateVector end = CoordinateVector.fromValues(1, 1, 1);
		final int polynomialDegree = 2;
		final QkQkSpace grid = new QkQkSpace(start, end,
		                                     Ints.asList(6, 6, 6));
		final TPVectorCellIntegral<ContinuousTPVectorFunction> valueValue =
			new TPVectorCellIntegral<>(TPVectorCellIntegral.ROT_ROT);
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
		final MixedRightHandSideIntegral<TPCell, ContinuousTPShapeFunction,
			ContinuousTPVectorFunction, QkQkFunction> rightHandSideIntegral =
			MixedRightHandSideIntegral.fromVelocityIntegral(
				new TPVectorRightHandSideIntegral<ContinuousTPVectorFunction>(MaxwellReferenceSolution.rightHandSide(),
				                                                              TPVectorRightHandSideIntegral.VALUE));
		final List<RightHandSideIntegral<TPCell, QkQkFunction>> rightHandSideIntegrals = new ArrayList<>();
		rightHandSideIntegrals.add(rightHandSideIntegral);
		
		final List<BoundaryRightHandSideIntegral<TPFace, QkQkFunction>> boundaryFaceIntegrals =
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
		final Vector solution1 = i.solveGMRES(new BlockSparseMatrix(grid.getSystemMatrix(), 20), grid.getRhs(),
		                                      1e-6);
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
		
		PlotWindow.addPlot(new MixedPlot3D(solut, grid.generatePlotPoints(20), 20));
		PlotWindow.addPlot(new MixedPlot3D(MaxwellReferenceSolution.mixedReferenceSolution(),
		                                   grid.generatePlotPoints(20),
		                                   20));
	}
}
