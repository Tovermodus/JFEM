package tensorproduct;

import basic.BoundaryRightHandSideIntegral;
import basic.Function;
import basic.ScalarShapeFunction;
import basic.VectorShapeFunction;
import linalg.CoordinateVector;

public class TPVectorBoundaryFaceIntegral<ST extends VectorShapeFunction<TPCell,TPFace,ST>> extends BoundaryRightHandSideIntegral<TPCell,TPFace,
	ST>
{
	static final String VALUE="Value";
	public TPVectorBoundaryFaceIntegral(Function<?,?,?> rightHandSide, String name, boolean weightIsTensorProduct)
	{
		super(rightHandSide, name);
		if(name.equals(VALUE) && !(rightHandSide.value(new CoordinateVector(rightHandSide.getDomainDimension())) instanceof Double))
			throw new IllegalArgumentException();
	}
	
	@Override
	public double evaluateBoundaryRightHandSideIntegral(TPFace face,
	                                                    ST shapeFunction1)
	{
		if(name.equals(VALUE))
		{
				return TPFaceIntegral.integrateNonTensorProduct(x->shapeFunction1.value(x).inner((CoordinateVector)(rightHandSide.value(x))),face.cell1Ds, face.flatDimension, face.otherCoordinate);
		}
		throw new UnsupportedOperationException("unknown rhs integral");
	}
	
}
