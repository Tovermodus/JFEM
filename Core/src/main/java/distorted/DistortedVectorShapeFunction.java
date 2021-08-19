package distorted;

import basic.SingleComponentVectorShapeFunction;
import distorted.geometry.DistortedCell;
import distorted.geometry.DistortedFace;
import linalg.CoordinateMatrix;
import linalg.CoordinateVector;

public class DistortedVectorShapeFunction extends SingleComponentVectorShapeFunction<DistortedCell, DistortedFace,
	DistortedShapeFunction, DistortedVectorShapeFunction>

{
	public DistortedVectorShapeFunction(final DistortedCell supportCell, final int polynomialDegree, final int localIndex)
	{
		super(supportCell, polynomialDegree, localIndex, DistortedShapeFunction.class);
	}
	
	public DistortedVectorShapeFunction(final DistortedShapeFunction componentFunction, final int component)
	{
		super(componentFunction, component);
	}
	
	public CoordinateVector valueOnReferenceCell(CoordinateVector x, DistortedCell cell)
	{
		return CoordinateVector.getUnitVector(getRangeDimension(), getComponent())
		                       .mul(getComponentFunction().valueOnReferenceCell(x, cell));
	}
	
	public CoordinateMatrix gradientOnReferenceCell(CoordinateVector x, DistortedCell cell)
	{
		return getComponentFunction().gradientOnReferenceCell(x, cell)
		                             .outer(CoordinateVector.getUnitVector(getDomainDimension(),
		                                                                   getComponent()));
	}
}
