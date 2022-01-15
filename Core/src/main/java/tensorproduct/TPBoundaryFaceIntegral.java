package tensorproduct;

import basic.BoundaryRightHandSideIntegral;
import basic.Function;
import basic.ScalarShapeFunction;
import linalg.CoordinateVector;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

public class TPBoundaryFaceIntegral<ST extends ScalarShapeFunction<TPCell, TPFace>>
	extends BoundaryRightHandSideIntegral<TPFace,
	ST>
{
	public static final String VALUE = "Value";
	public static final String VALUE_NORMAL = "ValueNormal";
	
	public TPBoundaryFaceIntegral(final Function<?, ?, ?> rightHandSide, final String name)
	{
		super(rightHandSide, name);
		if (name.equals(VALUE) && !(rightHandSide.value(new CoordinateVector(rightHandSide.getDomainDimension())) instanceof Double))
			throw new IllegalArgumentException();
		if (name.equals(VALUE_NORMAL) && !(rightHandSide.value(new CoordinateVector(rightHandSide.getDomainDimension())) instanceof Double))
			throw new IllegalArgumentException();
	}
	
	@Override
	public double evaluateBoundaryRightHandSideIntegral(final TPFace face,
	                                                    final ST shapeFunction1)
	{
		if (!face.isBoundaryFace())
			return 0;
		if (name.equals(VALUE))
		{
			return TPFaceIntegral.integrateNonTensorProduct(x -> shapeFunction1.value(x)
				                                                * (Double) (rightHandSide.value(x)),
			                                                face,
			                                                quadratureRule1D);
		}
		if (name.equals(VALUE_NORMAL))
		{
			return TPFaceIntegral.integrateNonTensorProduct(x -> (Double) rightHandSide.value(x)
				                                                * ((CoordinateVector) shapeFunction1.gradient(x)).inner(face.getNormal()
				                                                                                                            .value(x)),
			                                                face,
			                                                quadratureRule1D);
		}
		throw new UnsupportedOperationException("unknown rhs integral");
	}
}
