package distorted;

import basic.SingleComponentVectorShapeFunction;
import distorted.geometry.DistortedCell;
import distorted.geometry.DistortedFace;
import linalg.CoordinateVector;
import linalg.Rank1CoordinateMatrix;

public class DistortedVectorShapeFunction
	extends SingleComponentVectorShapeFunction<DistortedCell, DistortedFace, DistortedShapeFunction, DistortedVectorShapeFunction>
	implements DistortedVectorFunction

{
	public DistortedVectorShapeFunction(final DistortedCell supportCell,
	                                    final int polynomialDegree,
	                                    final int localIndex)
	{
		super(supportCell, polynomialDegree, localIndex, DistortedShapeFunction.class);
	}
	
	public DistortedVectorShapeFunction(final DistortedShapeFunction componentFunction, final int component)
	{
		super(componentFunction, component);
	}
	
	@Override
	public CoordinateVector valueOnReferenceCell(final CoordinateVector x, final DistortedCell cell)
	{
		return CoordinateVector
			.getUnitVector(getRangeDimension(), getComponent())
			.mul(getComponentFunction().valueOnReferenceCell(x, cell));
	}
	
	@Override
	public Rank1CoordinateMatrix gradientOnReferenceCell(final CoordinateVector x, final DistortedCell cell)
	{
		return getComponentFunction()
			.gradientOnReferenceCell(x, cell)
			.outer(CoordinateVector.getUnitVector(getDomainDimension(), getComponent()));
	}
}
