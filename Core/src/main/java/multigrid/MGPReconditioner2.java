package multigrid;

import basic.MetricWindow;
import linalg.DenseVector;
import linalg.IterativeSolverConvergenceMetric;
import linalg.Vector;
import linalg.VectorMultiplyable;

public class MGPReconditioner2
	implements VectorMultiplyable
{
	final MGInterface mg;
	IterativeSolverConvergenceMetric ocm;
	
	public MGPReconditioner2(final MGInterface mg)
	{
		this.mg = mg;
		ocm = new IterativeSolverConvergenceMetric(1e-2);
		MetricWindow.getInstance()
		            .setMetric("mgp", ocm);
	}
	
	@Override
	public int getVectorSize()
	{
		return mg.getFinestSystem()
		         .getVectorSize();
	}
	
	@Override
	public int getTVectorSize()
	{
		return mg.getFinestSystem()
		         .getVectorSize();
	}
	
	@Override
	public Vector mvMul(final Vector vector)
	{
		System.out.println();
		Vector v = new DenseVector(vector);
		final Vector residual = vector.sub(mg.getFinestSystem()
		                                     .mvMul(vector));
		final Vector residualSolution = mg.fullVCycleSolver(residual);
		v = vector.add(residualSolution);
		System.out.print("     mgp1  "
			                 + mg.getFinestSystem()
			                     .mvMul(vector)
			                     .sub(vector)
			                     .euclidianNorm());
		System.out.print("  mgp2  "
			                 + mg.getFinestSystem()
			                     .mvMul(v)
			                     .sub(vector)
			                     .euclidianNorm());
//		residual = vector.sub(mg.getFinestSystem()
//		                        .mvMul(v));
//		residualSolution = mg.fullVCycle(residual);
//		residualSolution = mg.vCycle(residualSolution, residual);
//		residualSolution = mg.vCycle(residualSolution, residual);
//		residualSolution = mg.vCycle(residualSolution, residual);
//		System.out.println("     mgp3     "
//			                   + mg.getFinestSystem()
//			                       .mvMul(v.add(residualSolution))
//			                       .sub(vector)
//			                       .euclidianNorm());
//		System.out.println();
		ocm.publishIterate(mg.getFinestSystem()
		                     .mvMul(vector)
		                     .sub(vector)
		                     .euclidianNorm() / mg.getFinestSystem()
		                                          .mvMul(v)
		                                          .sub(vector)
		                                          .euclidianNorm());
		return v;
	}
	
	@Override
	public Vector tvMul(final Vector vector)
	{
		throw new UnsupportedOperationException("not implemented yet");
	}
}
