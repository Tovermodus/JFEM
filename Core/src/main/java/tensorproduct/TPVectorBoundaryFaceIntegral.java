package tensorproduct;

import basic.BoundaryRightHandSideIntegral;
import basic.Function;
import basic.VectorShapeFunction;
import com.google.common.collect.Iterables;
import linalg.CoordinateVector;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

public class TPVectorBoundaryFaceIntegral<ST extends VectorShapeFunction<TPCell, TPFace>> extends BoundaryRightHandSideIntegral<TPFace,
	ST>
{
	public static final String VALUE="Value";
	public static final String NORMAL_VALUE = "NormalValue";
	
	
	public TPVectorBoundaryFaceIntegral(Function<?, ?, ?> rightHandSide, String name)
	{
		super(rightHandSide, name);
		if(name.equals(VALUE) && !(rightHandSide.value(new CoordinateVector(rightHandSide.getDomainDimension())) instanceof CoordinateVector))
			throw new IllegalArgumentException();
		if(name.equals(NORMAL_VALUE) && !(rightHandSide.value(new CoordinateVector(rightHandSide.getDomainDimension())) instanceof Double))
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
			return TPFaceIntegral.integrateNonTensorProduct(x->shapeFunction1.value(x).inner((CoordinateVector)(rightHandSide.value(x))),
				face, quadratureRule1D);
		}
		if(name.equals(NORMAL_VALUE))
		{
			TPCell cell = Iterables.getFirst(face.getCells(), null);
			if(cell == null)
				throw new IllegalStateException("face has no cells");
			return TPFaceIntegral.integrateNonTensorProduct(x->shapeFunction1.value(x).inner(cell.getOuterNormal(face).value(x))*(Double)rightHandSide.value(x),
				face, quadratureRule1D);
		}
		throw new UnsupportedOperationException("unknown rhs integral");
	}
	
}
