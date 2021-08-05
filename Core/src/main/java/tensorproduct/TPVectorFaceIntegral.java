package tensorproduct;

import basic.*;
import linalg.CoordinateVector;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

public class TPVectorFaceIntegral<ST extends VectorShapeFunction<TPCell, TPFace>> extends FaceIntegral<TPFace,ST>
{
	public static String VALUE_NORMALAVERAGE_GRAD_AVERAGE = "ValueNormalaverageGradAverage";
	public static String GRAD_AVERAGE_VALUE_NORMALAVERAGE = "GradAverageValueNormalaverage";
	public static String VALUE_NORMALAVERAGE_VALUE_NORMALAVERAGE = "ValueNormalaverageValueNormalaverage";
	public static String BOUNDARY_VALUE = "BoundaryValue";
	
	public TPVectorFaceIntegral(Function<?, ?, ?> weight, String name)
	{
		this(weight, name, QuadratureRule1D.Gauss5);
	}
	
	public TPVectorFaceIntegral(Function<?, ?, ?> weight, String name, QuadratureRule1D quadratureRule1D)
	{
		super(weight,name, quadratureRule1D);
		if(name.equals(VALUE_NORMALAVERAGE_GRAD_AVERAGE) && !(weight.value(new CoordinateVector(weight.getDomainDimension())) instanceof Double))
			throw new IllegalArgumentException();
		if(name.equals(GRAD_AVERAGE_VALUE_NORMALAVERAGE) && !(weight.value(new CoordinateVector(weight.getDomainDimension())) instanceof Double))
			throw new IllegalArgumentException();
		if(name.equals(VALUE_NORMALAVERAGE_VALUE_NORMALAVERAGE) && !(weight.value(new CoordinateVector(weight.getDomainDimension())) instanceof Double))
			throw new IllegalArgumentException();
	}
	public TPVectorFaceIntegral(String name, QuadratureRule1D quadratureRule1D)
	{
		super(name, quadratureRule1D);
	}
	public TPVectorFaceIntegral(String name)
	{
		this(name, QuadratureRule1D.Gauss5);
	}
	@Override
	public double evaluateFaceIntegral(TPFace face, ST shapeFunction1,
	                                   ST shapeFunction2)
	{
		if (name.equals(VALUE_NORMALAVERAGE_GRAD_AVERAGE))
		{
			return TPFaceIntegral.integrateNonTensorProduct(x -> shapeFunction1.normalAverageInValue(face
				,x).frobeniusInner(shapeFunction2.averageInDerivative(face,x))*(double)weight.value(x),
				face,
				quadratureRule1D);
		}
		if (name.equals(GRAD_AVERAGE_VALUE_NORMALAVERAGE))
		{
			return TPFaceIntegral.integrateNonTensorProduct(x -> shapeFunction2.normalAverageInValue(face
				,x).frobeniusInner(shapeFunction1.averageInDerivative(face,x))*(double)weight.value(x),
				face,
				quadratureRule1D);
		}
		if (name.equals(VALUE_NORMALAVERAGE_VALUE_NORMALAVERAGE))
		{
			return TPFaceIntegral.integrateNonTensorProduct(x -> shapeFunction1.normalAverageInValue(face
				,x).frobeniusInner(shapeFunction2.normalAverageInValue(face,x))*(double)weight.value(x),
				face,
				quadratureRule1D);
		}
		if (name.equals(BOUNDARY_VALUE))
		{
			if (face.isBoundaryFace())
				return TPFaceIntegral.integrateNonTensorProduct(x -> shapeFunction1.normalAverageInValue(face
					, x).frobeniusInner(shapeFunction2.normalAverageInValue(face, x)) * (double) weight.value(x),
					face,
					quadratureRule1D);
			else
				return 0;
		}
		throw new UnsupportedOperationException("unkown face integral name");
	}
	
}
