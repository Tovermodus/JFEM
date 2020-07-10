package tensorproduct;
import basic.*;
import linalg.CoordinateVector;

public class TPBoundaryFaceIntegral extends BoundaryRightHandSideIntegral<TPCell<TPShapeFunction>,TPFace<TPShapeFunction>,
	TPShapeFunction>
{
	static final String VALUE="Value";
	private final boolean weightIsTensorProduct;
	public TPBoundaryFaceIntegral(Function<?,?,?> rightHandSide, String name, boolean weightIsTensorProduct)
	{
		super(rightHandSide, name);
		if(name.equals(VALUE) && !(rightHandSide.value(new CoordinateVector(rightHandSide.getDomainDimension())) instanceof Double))
			throw new IllegalArgumentException();
		this.weightIsTensorProduct = weightIsTensorProduct;
	}
	
	@Override
	public double evaluateBoundaryRightHandSideIntegral(TPFace<TPShapeFunction> face,
	                                                    TPShapeFunction shapeFunction1)
	{
		if(name.equals(VALUE))
		{
			if(weightIsTensorProduct)
				return TPFaceIntegral.integrateTensorProduct(x->shapeFunction1.value(x)*(Double)(rightHandSide.value(x)),face.cell1Ds, face.flatDimension, face.otherCoordinate);
			else
				return TPFaceIntegral.integrateNonTensorProduct(x->shapeFunction1.value(x)*(Double)(rightHandSide.value(x)),face.cell1Ds, face.flatDimension, face.otherCoordinate);
		}
		throw new UnsupportedOperationException("unknown rhs integral");
	}
	
}
