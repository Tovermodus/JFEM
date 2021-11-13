package basic;

import linalg.AddableVectorMultiplyable;
import linalg.DenseMatrix;
import linalg.DenseVector;
import scala.Function2;

import java.util.ArrayList;
import java.util.function.Function;

public abstract class FullyImplicitHyperbolicIntegrator<M extends AddableVectorMultiplyable<M>>
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
		final M twoDerivOp = getDoubleDerivativeOperator();
		final M oneDerivOp = getSingleDerivativeOperator();
		final M zeroDerivOp = getNoDerivativeOperator();
		final M implicitOperator
			= twoDerivOp.mul(1. / (dt * dt))
			                .add(oneDerivOp.mul((1. / dt)))
			                .add(zeroDerivOp);
		DenseVector rhs = getAdditionalRhs();
		rhs.addInPlace(twoDerivOp.mvMul(currentIterate).mul(2.0/(dt*dt)));
		rhs.addInPlace(twoDerivOp.mvMul(lastIterate).mul(-1.0/(dt*dt)));
		rhs.addInPlace(oneDerivOp.mvMul(currentIterate).mul(1.0/(dt)));
		lastIterate = currentIterate;
		currentIterate = getSolver().apply(implicitOperator, rhs);
	}
	
	protected abstract DenseVector getAdditionalRhs();
	
	protected abstract M getDoubleDerivativeOperator();
	
	protected abstract M getSingleDerivativeOperator();
	
	protected abstract M getNoDerivativeOperator();
	
	protected abstract Function2<M, DenseVector, DenseVector> getSolver();
	
	protected abstract void postIterationCallback();
}
