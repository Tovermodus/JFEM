package basic;

import linalg.AddableVectorMultiplyable;
import linalg.DenseMatrix;
import linalg.DenseVector;
import scala.Function2;

public abstract class FullyImplicitParabolicIntegrator<M extends AddableVectorMultiplyable<M>>
{
	final public double dt;
	final public int timeSteps;
	protected DenseVector currentIterate;
	final public DenseMatrix iterateHistory;
	protected double time = 0.0;
	
	public FullyImplicitParabolicIntegrator(final double dt, final int timeSteps)
	{
		this.dt = dt;
		this.timeSteps = timeSteps;
		currentIterate = initializeInitialIterate();
		iterateHistory = new DenseMatrix(timeSteps, currentIterate.getLength());
		iterateHistory.addRow(currentIterate, 0);
	}
	
	protected abstract DenseVector initializeInitialIterate();
	
	public void loop()
	{
		for (int i = 0; i < timeSteps; i++)
		{
			time += dt;
			timeStep();
			postIterationCallback();
		}
	}
	
	private void timeStep()
	{
		final M oneDerivOp = getSingleDerivativeOperator();
		final M zeroDerivOp = getNoDerivativeOperator();
		final M implicitOperator
			= oneDerivOp.mul((1. / dt))
			            .add(zeroDerivOp);
		final DenseVector rhs = getAdditionalRhs();
		rhs.addInPlace(oneDerivOp.mvMul(currentIterate)
		                         .mul(1.0 / (dt)));
		currentIterate = getSolver().apply(implicitOperator, rhs);
	}
	
	protected abstract DenseVector getAdditionalRhs();
	
	protected abstract M getSingleDerivativeOperator();
	
	protected abstract M getNoDerivativeOperator();
	
	protected abstract Function2<M, DenseVector, DenseVector> getSolver();
	
	protected abstract void postIterationCallback();
}
