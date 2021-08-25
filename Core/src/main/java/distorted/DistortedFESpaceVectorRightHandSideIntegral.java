package distorted;

import basic.RightHandSideIntegral;
import distorted.geometry.DistortedCell;

public class DistortedFESpaceVectorRightHandSideIntegral extends RightHandSideIntegral<DistortedCell, DistortedVectorShapeFunction>
{
	public static final String H1 = "H1";
	public static final String VALUE = "Value";
	DistortedVectorFunction rightHandSide;
	
	public DistortedFESpaceVectorRightHandSideIntegral(final DistortedVectorFunction rightHandSide, final String name)
	{
		super(rightHandSide, name);
		this.rightHandSide = rightHandSide;
	}
	
	@Override
	public double evaluateRightHandSideIntegral(final DistortedCell cell, final DistortedVectorShapeFunction shapeFunction1)
	{
		if (name.equals(H1))
		{
			return DistortedCellIntegral.integrateOnReferenceCell(x -> shapeFunction1
				                                                      .valueOnReferenceCell(x, cell)
				                                                      .inner(rightHandSide.valueOnReferenceCell(x, cell)) + shapeFunction1
				                                                      .gradientOnReferenceCell(x, cell)
				                                                      .frobeniusInner(rightHandSide.gradientOnReferenceCell(x, cell)), cell,
			                                                      quadratureRule1D);
		}
		if (name.equals(VALUE))
		{
			return DistortedCellIntegral.integrateOnReferenceCell(x -> shapeFunction1
				.valueOnReferenceCell(x, cell)
				.inner(rightHandSide.valueOnReferenceCell(x, cell)), cell, quadratureRule1D);
		} else throw new IllegalArgumentException("Name unknown");
	}
}
