package dlm;

import basic.VectorFunctionOnCells;
import com.google.common.base.Stopwatch;
import it.unimi.dsi.fastutil.ints.Int2DoubleMap;
import linalg.*;
import mixed.TaylorHoodSpace;
import org.jetbrains.annotations.NotNull;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

import java.util.List;

public abstract class MultiGridFluid
	implements Fluid
{
	final TaylorHoodSpace space;
	final CoordinateVector startCoordinates;
	final CoordinateVector endCoordinates;
	final IntCoordinates coarsestCells;
	final int refinements;
	final int polynomialDegree;
	
	public MultiGridFluid(final CoordinateVector startCoordinates, final CoordinateVector endCoordinates,
	                      final IntCoordinates coarsestCells, final int refinements, final int polynomialDegree)
	{
		this.startCoordinates = startCoordinates;
		this.endCoordinates = endCoordinates;
		this.coarsestCells = coarsestCells;
		this.refinements = refinements;
		this.polynomialDegree = polynomialDegree;
		space = new TaylorHoodSpace(startCoordinates,
		                            endCoordinates,
		                            coarsestCells.mul((int) Math.pow(2, refinements)));
		space.assembleCells();
		space.assembleFunctions(polynomialDegree);
	}
	
	@Override
	public TaylorHoodSpace getSpace()
	{
		return space;
	}
	
	@Override
	public FluidSystem buildSystem(final double t, final FluidIterate iterate, final List<Particle> particles)
	{
		return getFluidSystemForSpace(getSpace(), getVelocity(iterate), t, iterate.current);
	}
	
	@NotNull
	FluidSystem getFluidSystemForSpace(final TaylorHoodSpace space,
	                                   final VectorFunctionOnCells<TPCell, TPFace> velocity,
	                                   final double t,
	                                   final Vector currentIterate)
	{
		final int n = space.getShapeFunctionMap()
		                   .size();
		final SparseMatrix massMatrix = new SparseMatrix(n, n);
		space.writeCellIntegralsToMatrix(getMassIntegrals(), massMatrix);
		final SparseMatrix flowMatrix = new SparseMatrix(n, n);
		space.writeCellIntegralsToMatrix(getIntegrals(), flowMatrix);
		final SparseMatrix semiImplicitMatrix = new SparseMatrix(n, n);
		System.out.println("semiImp");
		final Stopwatch s = Stopwatch.createStarted();
		space.writeCellIntegralsToMatrix(getSemiImplicitIntegrals(velocity),
		                                 semiImplicitMatrix);
		System.out.println("semiImpDone" + s.elapsed());
		final DenseVector forceRhs = new DenseVector(n);
		space.writeCellIntegralsToRhs(getForceIntegrals(t), forceRhs);
		final DenseVector accelRhs = massMatrix.mvMul(currentIterate);
		return new FluidSystem(massMatrix, flowMatrix, semiImplicitMatrix, forceRhs, accelRhs);
	}
}
