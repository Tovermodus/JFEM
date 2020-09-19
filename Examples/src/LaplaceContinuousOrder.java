import basic.*;
import com.google.common.base.Stopwatch;
import com.google.common.primitives.Ints;
import linalg.CoordinateVector;
import linalg.IterativeSolver;
import linalg.SparseMatrix;
import linalg.Vector;
import tensorproduct.*;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.TreeSet;

public class LaplaceContinuousOrder
{
	public static void main(String[] args)
	{
		
		CoordinateVector start = CoordinateVector.fromValues(-1, -1);
		CoordinateVector end = CoordinateVector.fromValues(1, 1);
		
		TPVectorCellIntegral<ContinuousTPVectorFunction> gg =
			new TPVectorCellIntegral<>(TPVectorCellIntegral.GRAD_GRAD);
		ArrayList<CellIntegral<TPCell, TPFace, ContinuousTPVectorFunction>> cellIntegrals =
			new ArrayList<>();
		cellIntegrals.add(gg);
		ArrayList<FaceIntegral<TPCell, TPFace, ContinuousTPVectorFunction>> faceIntegrals = new ArrayList<>();
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
					return CoordinateVector.fromValues(0,
						0);
				}
			},
				TPVectorRightHandSideIntegral.VALUE);
		ArrayList<RightHandSideIntegral<TPCell, TPFace, ContinuousTPVectorFunction>> rightHandSideIntegrals = new ArrayList<>();
		rightHandSideIntegrals.add(rightHandSideIntegral);
		ArrayList<BoundaryRightHandSideIntegral<TPCell, TPFace, ContinuousTPVectorFunction>> boundaryFaceIntegrals = new ArrayList<>();
		
		int polynomialDegree = 3;
		List<ScalarFunction> solutions = new ArrayList<>();
		List<VectorFunction> solutionsVec = new ArrayList<>();
		ContinuousTPFEVectorSpace grid = null;
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
					return 2*(1+pos.y())/((3+pos.x())*(3+pos.x())+(1+pos.y())*(1+pos.y()));
			}
		};
		for(int i = 0; i < 5; i++)
		{
			grid = new ContinuousTPFEVectorSpace(start, end,
				Ints.asList(2*(int)Math.pow(2,i),2*(int)Math.pow(2,i)), polynomialDegree);
			grid.assembleCells();
			grid.assembleFunctions(polynomialDegree);
			grid.initializeSystemMatrix();
			grid.initializeRhs();
			grid.evaluateCellIntegrals(cellIntegrals, rightHandSideIntegrals);
			grid.evaluateFaceIntegrals(faceIntegrals, boundaryFaceIntegrals);
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
					if (Math.abs(pos.x()) == 1||Math.abs(pos.y()) == 1)
						return CoordinateVector.fromValues(func.value(pos), func.value(pos));
					return CoordinateVector.fromValues(0,0);
				}
			});
			IterativeSolver<SparseMatrix> it = new IterativeSolver<>();
			Vector solution1 = it.solveGMRES(grid.getSystemMatrix(), grid.getRhs(), 1e-9);
			VectorFESpaceFunction<ContinuousTPVectorFunction> solut =
				new VectorFESpaceFunction<>(
					grid.getShapeFunctions(), solution1);
			solutionsVec.add(new VectorFunction()
			{
				@Override
				public int getDomainDimension()
				{
					return solut.getDomainDimension();
				}
				
				@Override
				public CoordinateVector value(CoordinateVector pos)
				{
					return solut.value(pos);
				}
			});
		}
		solutionsVec.add(new VectorFunction()
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
		System.out.println(ConvergenceOrderEstimator.estimateL2Vector(solutionsVec,
			grid.generatePlotPoints(40)));
		
	}
}
