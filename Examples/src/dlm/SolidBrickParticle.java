package dlm;

import basic.CellIntegral;
import basic.RightHandSideIntegral;
import distorted.*;
import distorted.geometry.DistortedCell;
import it.unimi.dsi.fastutil.ints.Int2DoubleArrayMap;
import it.unimi.dsi.fastutil.ints.Int2DoubleMap;
import linalg.CoordinateMatrix;
import linalg.CoordinateTensor;
import linalg.CoordinateVector;

import java.util.ArrayList;
import java.util.List;
import java.util.function.Function;

public class SolidBrickParticle
	implements Particle
{
	private final CoordinateVector initialVelocity;
	private final double lameLambda;
	private final double lameMu;
	BrickVectorSpace space;
	private final double density;
	
	public SolidBrickParticle(final CoordinateVector center,
	                          final double width,
	                          final double height,
	                          final int refines,
	                          final int polynomialDegree,
	                          final CoordinateVector initialVelocity,
	                          final double lameLambda,
	                          final double lameMu, final double density)
	{
		this.initialVelocity = initialVelocity;
		this.lameLambda = lameLambda;
		this.lameMu = lameMu;
		this.density = density;
		space = new BrickVectorSpace(center, width, height, refines);
		space.assembleCells();
		space.assembleFunctions(polynomialDegree);
	}
	
	@Override
	public DistortedGridSpace<DistortedVectorShapeFunction, CoordinateVector, CoordinateMatrix, CoordinateTensor> getSpace()
	{
		return space;
	}
	
	@Override
	public List<CellIntegral<DistortedCell, DistortedVectorShapeFunction>> getIntegrals()
	{
		final DistortedVectorCellIntegral lambdaIntegral = new DistortedVectorCellIntegral(lameLambda,
		                                                                                   DistortedVectorCellIntegral.TRACE_SYM);
		final DistortedVectorCellIntegral muIntegral = new DistortedVectorCellIntegral(lameMu * 2,
		                                                                               DistortedVectorCellIntegral.SYM_SYM);
		return List.of(lambdaIntegral, muIntegral);
	}
	
	@Override
	public List<CellIntegral<DistortedCell, DistortedVectorShapeFunction>> getLagrangeIntegrals()
	{
		final DistortedVectorCellIntegral l2Integral = new DistortedVectorCellIntegral(1,
		                                                                               DistortedVectorCellIntegral.VALUE_VALUE);
		return List.of(l2Integral);
	}
	
	@Override
	public List<CellIntegral<DistortedCell, DistortedVectorShapeFunction>> getMassIntegrals()
	{
		final DistortedVectorCellIntegral massIntegral = new DistortedVectorCellIntegral(density,
		                                                                                 DistortedVectorCellIntegral.VALUE_VALUE);
		return List.of(massIntegral);
	}
	
	@Override
	public List<RightHandSideIntegral<DistortedCell, DistortedVectorShapeFunction>> getForceIntegrals(final double t)
	{
		return new ArrayList<>();
	}
	
	@Override
	public List<RightHandSideIntegral<DistortedCell, DistortedVectorShapeFunction>> getBackgroundLagrangeIntegrals(
		final DistortedVectorFunctionOnCells backgroundFunctionAtX)
	{
		final DistortedVectorDistortedRightHandSideIntegral lagrange =
			new DistortedVectorDistortedRightHandSideIntegral(backgroundFunctionAtX,
			                                                  DistortedVectorDistortedRightHandSideIntegral.VALUE);
		return List.of(lagrange);
	}
	
	@Override
	public List<CellIntegral<DistortedCell, DistortedVectorShapeFunction>> getSemiImplicitIntegrals(
		final DistortedVectorFunctionOnCells displacement)
	{
		return new ArrayList<>();
	}
	
	@Override
	public Function<CoordinateVector, CoordinateVector> getInitialVelocity()
	{
		return x -> initialVelocity;
	}
	
	@Override
	public Int2DoubleMap getDirichletNodeValues(final double t)
	{
		return new Int2DoubleArrayMap();
	}
}
