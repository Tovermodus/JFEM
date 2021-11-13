package basic;

import linalg.DenseMatrix;
import linalg.DenseVector;
import linalg.VectorMultiplyable;
import scala.Function2;

public abstract class FullyImplicitParabolicIntegrator
{
	final public double dt;
	final public int timeSteps;
	protected DenseVector currentIterate;
	public DenseMatrix iterateHistory;
	protected double time = 0.0;
	
	public FullyImplicitParabolicIntegrator(final double dt, final int timeSteps)
	{
		this.dt = dt;
		this.timeSteps = timeSteps;
	}
	
	protected abstract DenseVector initializeInitialIterate();
	
	public void loop()
	{
		currentIterate = initializeInitialIterate();
		iterateHistory = new DenseMatrix(timeSteps + 2, currentIterate.getLength());
		iterateHistory.addRow(currentIterate, 0);
		for (int i = 0; i < timeSteps; i++)
		{
			time += dt;
			timeStep();
			iterateHistory.addRow(currentIterate, i + 2);
			postIterationCallback();
		}
	}
	
	private void timeStep()
	{
		final VectorMultiplyable oneDerivOp = getSingleDerivativeOperator();
		final VectorMultiplyable zeroDerivOp = getNoDerivativeOperator();
		final VectorMultiplyable implicitOperator
			= oneDerivOp.mul((1. / dt))
			            .addVm(zeroDerivOp);
		final DenseVector rhs = getAdditionalRhs();
		rhs.addInPlace(oneDerivOp.mvMul(currentIterate)
		                         .mul(1.0 / (dt)));
		currentIterate = getSolver().apply(implicitOperator, rhs);
	}
	
	protected abstract DenseVector getAdditionalRhs();
	
	protected abstract VectorMultiplyable getSingleDerivativeOperator();
	
	protected abstract VectorMultiplyable getNoDerivativeOperator();
	
	protected abstract Function2<VectorMultiplyable, DenseVector, DenseVector> getSolver();
	
	protected abstract void postIterationCallback();
}
