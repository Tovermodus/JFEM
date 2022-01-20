package dlm;

import com.google.common.base.Stopwatch;
import linalg.CoordinateVector;
import linalg.DenseVector;
import linalg.IntCoordinates;
import linalg.SparseMatrix;
import mixed.TaylorHoodSpace;

public abstract class SingleGridFluid
	implements Fluid
{
	private final TaylorHoodSpace space;
	
	public SingleGridFluid(final CoordinateVector startCoordinates, final CoordinateVector endCoordinates,
	                       final IntCoordinates cells, final int polynomialDegree)
	{
		space = new TaylorHoodSpace(startCoordinates, endCoordinates, cells);
		space.assembleCells();
		space.assembleFunctions(polynomialDegree);
	}
	
	@Override
	public TaylorHoodSpace getSpace()
	{
		return space;
	}
	
	@Override
	public FluidSystem buildSystem(final double t, final FluidIterate iterate)
	{
		final SparseMatrix massMatrix = new SparseMatrix(getSystemSize(), getSystemSize());
		getSpace().writeCellIntegralsToMatrix(getMassIntegrals(), massMatrix);
		final SparseMatrix flowMatrix = new SparseMatrix(getSystemSize(), getSystemSize());
		getSpace().writeCellIntegralsToMatrix(getIntegrals(), flowMatrix);
		final SparseMatrix semiImplicitMatrix = new SparseMatrix(getSystemSize(), getSystemSize());
		System.out.println("semiImp");
		final Stopwatch s = Stopwatch.createStarted();
		getSpace().writeCellIntegralsToMatrix(getSemiImplicitIntegrals(getVelocity(iterate)),
		                                      semiImplicitMatrix);
		System.out.println("semiImpDone" + s.elapsed());
		final DenseVector forceRhs = new DenseVector(getSystemSize());
		getSpace().writeCellIntegralsToRhs(getForceIntegrals(t), forceRhs);
		final DenseVector accelRhs = massMatrix.mvMul(iterate.current);
		return new FluidSystem(massMatrix, flowMatrix, semiImplicitMatrix, forceRhs, accelRhs);
	}
}
