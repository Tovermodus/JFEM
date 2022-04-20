package dlm;

import basic.CellIntegral;
import basic.LagrangeNodeFunctional;
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

public class FixedSolidParticle
	implements Particle
{
	private final double lameLambda;
	private final double lameMu;
	CircleVectorSpace space;
	private final double density;
	
	public FixedSolidParticle(final CoordinateVector center,
	                          final double radius,
	                          final int refines,
	                          final int polynomialDegree,
	                          final double lameLambda,
	                          final double lameMu, final double density)
	{
		this.lameLambda = lameLambda;
		this.lameMu = lameMu;
		this.density = density;
		space = new CircleVectorSpace(center, radius, refines);
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
		                                                                               DistortedVectorCellIntegral.H1);
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
			                                                  DistortedVectorDistortedRightHandSideIntegral.H1);
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
		return x -> new CoordinateVector(2);
	}
	
	@Override
	public Int2DoubleMap getDirichletNodeValues(final double t)
	{
		final Int2DoubleMap ret = new Int2DoubleArrayMap();
		for (final DistortedVectorShapeFunction sf : space.getShapeFunctions())
		{
			final CoordinateVector functionalPoint = ((LagrangeNodeFunctional) sf.getNodeFunctional()
			                                                                     .getComponentNodeFunctional()).getPoint();
			if (functionalPoint.dist(space.center) < space.radius / 3)
			{
				ret.put(sf.getGlobalIndex(), 0);
			}
		}
		return ret;
	}
}
