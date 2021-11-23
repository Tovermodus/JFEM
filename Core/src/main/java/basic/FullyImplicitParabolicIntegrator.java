package basic;

import linalg.DenseMatrix;
import linalg.DenseVector;
import linalg.Matrix;
import linalg.MutableVector;
import scala.Function2;

import java.util.function.BiConsumer;

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
		final Matrix oneDerivOp = getSingleDerivativeOperator();
		final Matrix zeroDerivOp = getNoDerivativeOperator();
		final Matrix implicitOperator
			= oneDerivOp.mul((1. / dt))
			            .add(zeroDerivOp);
		final DenseVector rhs = getAdditionalRhs();
		rhs.addInPlace(oneDerivOp.mvMul(currentIterate)
		                         .mul(1.0 / (dt)));
		boundaryApplier().accept(implicitOperator, rhs);
		currentIterate = getSolver().apply(implicitOperator, rhs);
	}
	
	protected abstract BiConsumer<Matrix, MutableVector> boundaryApplier();
	
	protected abstract DenseVector getAdditionalRhs();
	
	protected abstract Matrix getSingleDerivativeOperator();
	
	protected abstract Matrix getNoDerivativeOperator();
	
	protected abstract Function2<Matrix, DenseVector, DenseVector> getSolver();
	
	protected abstract void postIterationCallback();
}
