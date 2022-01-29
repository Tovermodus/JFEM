package dlm;

import linalg.*;
import mixed.*;
import multigrid.MGPreconditionerSpace;
import scala.Function2;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

public class TestIterativeImplicitSchur
	extends ImplicitSchurSolver
{
	final IterativeSolver it = new IterativeSolver(true);
	private final MGPreconditionerSpace<TaylorHoodSpace, TPCell, TPFace, QkQkFunction, MixedValue, MixedGradient
		, MixedHessian> mg;
	
	public TestIterativeImplicitSchur(final BlockSparseMatrix blockMatrix,
	                                  final MGPreconditionerSpace<TaylorHoodSpace, TPCell, TPFace, QkQkFunction, MixedValue, MixedGradient
		                                  , MixedHessian> mg)
	{
		super(blockMatrix);
		this.mg = mg;
		//firstInverse = ((SparseMatrix) mg.getFinestSystem()).inverse();
	}
	
	@Override
	protected Function2<VectorMultiplyable, Vector, Vector> solveSchur()
	{
		return (A, b) ->
		{
			it.showProgress = true;
			it.gm.ITERATIONS_BEFORE_RESTART = 30;
			it.gm.MAX_RESTARTS = 2;
			final Vector sol1 = it.solvePGMRES(A, mg, b, 1e-10);
			//final Vector sol2 = it2.solvePGMRES(A, firstInverse, b, 1e-10);
			System.out.println("ITERATIONS " + it.iterations);
			it.iterations = 0;
			System.out.println("sol1 done");
			System.out.println("res " + A.mvMul(sol1)
			                             .sub(b)
			                             .euclidianNorm());
			return sol1;
		};
	}
}
