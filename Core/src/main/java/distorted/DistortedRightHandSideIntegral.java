package distorted;

import basic.RightHandSideIntegral;
import basic.ScalarFunction;
import distorted.geometry.DistortedCell;

public class DistortedRightHandSideIntegral extends RightHandSideIntegral<DistortedCell, DistortedShapeFunction>
{
	public static final String VALUE = "VALUE";
	
	public DistortedRightHandSideIntegral(final ScalarFunction rightHandSide, final String name)
	{
		super(rightHandSide, name);
	}
	
	@Override
	public double evaluateRightHandSideIntegral(final DistortedCell cell, final DistortedShapeFunction shapeFunction1)
	{
		if (name.equals(VALUE))
		{
			return DistortedCellIntegral.integrateOnReferenceCell(x -> shapeFunction1
				                                                      .valueOnReferenceCell(x, cell)
				                                                      * (Double) rightHandSide.value(cell.transformFromReferenceCell(x)),
			                                                      cell,
			                                                      quadratureRule1D);
		} else throw new IllegalArgumentException("Name unknown");
	}
}
