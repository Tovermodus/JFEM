package tensorproduct;
import basic.*;
import linalg.CoordinateVector;

public class TPBoundaryFaceIntegral<ST extends ScalarShapeFunction<TPCell,TPFace,TPEdge>> extends BoundaryRightHandSideIntegral<TPFace,
	ST>
{
	public static final String VALUE="Value";
	public static final String VALUE_NORMAL = "ValueNormal";
	public TPBoundaryFaceIntegral(Function<?,?,?> rightHandSide, String name)
	{
		super(rightHandSide, name);
		if(name.equals(VALUE) && !(rightHandSide.value(new CoordinateVector(rightHandSide.getDomainDimension())) instanceof Double))
			throw new IllegalArgumentException();
		if(name.equals(VALUE_NORMAL) && !(rightHandSide.value(new CoordinateVector(rightHandSide.getDomainDimension())) instanceof CoordinateVector))
			throw new IllegalArgumentException();
	}
	
	@Override
	public double evaluateBoundaryRightHandSideIntegral(TPFace face,
	                                                    ST shapeFunction1)
	{
		if(!face.isBoundaryFace())
			return 0;
		if(name.equals(VALUE))
		{
			return TPFaceIntegral.integrateNonTensorProduct(x->shapeFunction1.value(x)
					*(Double)(rightHandSide.value(x)),
				face.cell1Ds,
				face.flatDimension,
				face.otherCoordinate,
				quadratureRule1D);
		}
		if(name.equals(VALUE_NORMAL))
		{
			return TPFaceIntegral.integrateNonTensorProduct(x->shapeFunction1.value(x)
				*(Double)(rightHandSide.value(x)),
				face.cell1Ds,
				face.flatDimension,
				face.otherCoordinate,
				quadratureRule1D);
		}
		throw new UnsupportedOperationException("unknown rhs integral");
	}
	
}
