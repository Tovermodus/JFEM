package tensorproduct;
import basic.*;
import linalg.CoordinateVector;

public class TPBoundaryFaceIntegral<ST extends ScalarShapeFunction<TPCell,TPFace,ST>> extends BoundaryRightHandSideIntegral<TPCell,TPFace,
	ST>
{
	public static final String VALUE="Value";
	private final boolean weightIsTensorProduct;
	public TPBoundaryFaceIntegral(Function<?,?,?> rightHandSide, String name, boolean weightIsTensorProduct)
	{
		super(rightHandSide, name);
		if(name.equals(VALUE) && !(rightHandSide.value(new CoordinateVector(rightHandSide.getDomainDimension())) instanceof Double))
			throw new IllegalArgumentException();
		this.weightIsTensorProduct = weightIsTensorProduct;
	}
	
	@Override
	public double evaluateBoundaryRightHandSideIntegral(TPFace face,
	                                                    ST shapeFunction1)
	{
		if(!face.isBoundaryFace())
			return 0;
		if(name.equals(VALUE))
		{
				return TPFaceIntegral.integrateNonTensorProduct(x->shapeFunction1.value(x)*(Double)(rightHandSide.value(x)),face.cell1Ds, face.flatDimension, face.otherCoordinate);
		}
		throw new UnsupportedOperationException("unknown rhs integral");
	}
	
}
