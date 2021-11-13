package basic;

import linalg.DenseMatrix;
import linalg.DenseVector;
import linalg.VectorMultiplyable;
import scala.Function2;

public abstract class FullyImplicitHyperbolicIntegrator
{
	final public double dt;
	final public int timeSteps;
	protected DenseVector currentIterate;
	protected DenseVector lastIterate;
	final public DenseMatrix iterateHistory;
	protected double time = 0.0;
	
	public FullyImplicitHyperbolicIntegrator(final double dt, final int timeSteps)
	{
		this.dt = dt;
		this.timeSteps = timeSteps;
		currentIterate = initializeInitialIterate();
		lastIterate = currentIterate.sub(initializeInitialDerivative().mul(dt));
		iterateHistory = new DenseMatrix(timeSteps, currentIterate.getLength());
		iterateHistory.addRow(lastIterate, 0);
		iterateHistory.addRow(currentIterate, 1);
	}
	
	protected abstract DenseVector initializeInitialIterate();
	
	protected abstract DenseVector initializeInitialDerivative();
	
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
		final VectorMultiplyable twoDerivOp = getDoubleDerivativeOperator();
		final VectorMultiplyable oneDerivOp = getSingleDerivativeOperator();
		final VectorMultiplyable zeroDerivOp = getNoDerivativeOperator();
		final VectorMultiplyable implicitOperator
			= twoDerivOp.mul(1. / (dt * dt))
			            .addVm(oneDerivOp.mul((1. / dt)))
			            .addVm(zeroDerivOp);
		final DenseVector rhs = getAdditionalRhs();
		rhs.addInPlace(twoDerivOp.mvMul(currentIterate)
		                         .mul(2.0 / (dt * dt)));
		rhs.addInPlace(twoDerivOp.mvMul(lastIterate)
		                         .mul(-1.0 / (dt * dt)));
		rhs.addInPlace(oneDerivOp.mvMul(currentIterate)
		                         .mul(1.0 / (dt)));
		lastIterate = currentIterate;
		currentIterate = getSolver().apply(implicitOperator, rhs);
	}
	
	protected abstract DenseVector getAdditionalRhs();
	
	protected abstract VectorMultiplyable getDoubleDerivativeOperator();
	
	protected abstract VectorMultiplyable getSingleDerivativeOperator();
	
	protected abstract VectorMultiplyable getNoDerivativeOperator();
	
	protected abstract Function2<VectorMultiplyable, DenseVector, DenseVector> getSolver();
	
	protected abstract void postIterationCallback();
}
