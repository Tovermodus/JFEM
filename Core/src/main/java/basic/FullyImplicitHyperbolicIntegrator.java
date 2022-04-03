package basic;

import io.vavr.Function2;
import linalg.*;

public abstract class FullyImplicitHyperbolicIntegrator
{
	final public double dt;
	final public int timeSteps;
	private MutableVector currentIterate;
	private MutableVector lastIterate;
	private DenseMatrix iterateHistory;
	private double time = 0.0;
	
	public FullyImplicitHyperbolicIntegrator(final double dt, final int timeSteps)
	{
		this.dt = dt;
		this.timeSteps = timeSteps;
	}
	
	protected abstract DenseVector initializeInitialIterate();
	
	protected abstract DenseVector initializeInitialDerivative();
	
	public void loop()
	{
		currentIterate = initializeInitialIterate();
		lastIterate = new DenseVector(currentIterate.sub(initializeInitialDerivative().mul(dt)));
		iterateHistory = new DenseMatrix(timeSteps + 2, currentIterate.getLength());
		iterateHistory.addRow(lastIterate, 0);
		iterateHistory.addRow(currentIterate, 1);
		postIterationCallback();
		for (int i = 0; i < timeSteps; i++)
		{
			time += dt;
			timeStep();
			postIterationCallback();
			iterateHistory.addRow(currentIterate, i + 2);
		}
	}
	
	private void timeStep()
	{
		final Matrix twoDerivOp = getDoubleDerivativeOperator();
		final Matrix oneDerivOp = getSingleDerivativeOperator();
		final Matrix zeroDerivOp = getNoDerivativeOperator();
		System.out.println("created Operators");
		Matrix implicitOperator
			= getOperator(twoDerivOp, oneDerivOp, zeroDerivOp);
		final DenseVector rhs = getAdditionalRhs();
		rhs.addInPlace(twoDerivOp.mvMul(currentIterate)
		                         .mul(2.0 / (dt * dt)));
		rhs.addInPlace(twoDerivOp.mvMul(lastIterate)
		                         .mul(-1.0 / (dt * dt)));
		rhs.addInPlace(oneDerivOp.mvMul(currentIterate)
		                         .mul(1.0 / (dt)));
		lastIterate = currentIterate;
		System.out.println("created all");
		implicitOperator = boundaryApplier().apply(implicitOperator, rhs);
		System.out.println("applied bdr");
		currentIterate = getSolver().apply(implicitOperator, rhs);
		System.out.println("solved");
	}
	
	public Matrix getOperator(final Matrix twoDerivOp, final Matrix oneDerivOp, final Matrix zeroDerivOp)
	{
		return twoDerivOp.mul(1. / (dt * dt))
		                 .add(oneDerivOp.mul((1. / dt)))
		                 .add(zeroDerivOp);
	}
	
	protected abstract Function2<Matrix, MutableVector, Matrix> boundaryApplier();
	
	protected abstract DenseVector getAdditionalRhs();
	
	protected abstract Matrix getDoubleDerivativeOperator();
	
	protected abstract Matrix getSingleDerivativeOperator();
	
	protected abstract Matrix getNoDerivativeOperator();
	
	protected abstract Function2<Matrix, Vector, MutableVector> getSolver();
	
	protected abstract void postIterationCallback();
	
	public Vector getCurrentIterate()
	{
		return currentIterate;
	}
	
	public Vector getLastIterate()
	{
		return lastIterate;
	}
	
	public Matrix getIterateHistory()
	{
		return iterateHistory;
	}
	
	public double getTime()
	{
		return time;
	}
}
