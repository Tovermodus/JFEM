package tensorproduct;

import basic.*;
import com.google.common.collect.BoundType;
import linalg.CoordinateVector;

import java.util.List;
import java.util.function.ToDoubleFunction;

public class TPVectorFaceIntegral<ST extends VectorShapeFunction<TPCell,TPFace,ST>> extends FaceIntegral<TPCell,
	TPFace,ST>
{
	public static String VALUE_NORMALAVERAGE_GRAD_AVERAGE = "ValueNormalaverageGradAverage";
	public static String GRAD_AVERAGE_VALUE_NORMALAVERAGE = "GradAverageValueNormalaverage";
	public static String VALUE_NORMALAVERAGE_VALUE_NORMALAVERAGE = "ValueNormalaverageValueNormalaverage";
	public static String BOUNDARY_VALUE = "BoundaryValue";
	public TPVectorFaceIntegral(Function<?,?,?> weight, String name)
	{
		super(weight,name);
		if(name.equals(VALUE_NORMALAVERAGE_GRAD_AVERAGE) && !(weight.value(new CoordinateVector(weight.getDomainDimension())) instanceof Double))
			throw new IllegalArgumentException();
		if(name.equals(GRAD_AVERAGE_VALUE_NORMALAVERAGE) && !(weight.value(new CoordinateVector(weight.getDomainDimension())) instanceof Double))
			throw new IllegalArgumentException();
		if(name.equals(VALUE_NORMALAVERAGE_VALUE_NORMALAVERAGE) && !(weight.value(new CoordinateVector(weight.getDomainDimension())) instanceof Double))
			throw new IllegalArgumentException();
	}
	public TPVectorFaceIntegral(String name)
	{
		super(name);
	}
	@Override
	public double evaluateFaceIntegral(TPFace face, ST shapeFunction1,
	                                   ST shapeFunction2)
	{
		if (name.equals(VALUE_NORMALAVERAGE_GRAD_AVERAGE))
		{
			return TPFaceIntegral.integrateNonTensorProduct(x -> shapeFunction1.normalAverageInValue(face
				,x).frobeniusInner(shapeFunction2.averageInDerivative(face,x))*(double)weight.value(x),
				face.cell1Ds,
				face.flatDimension,
				face.otherCoordinate);
		}
		if (name.equals(GRAD_AVERAGE_VALUE_NORMALAVERAGE))
		{
			return TPFaceIntegral.integrateNonTensorProduct(x -> shapeFunction2.normalAverageInValue(face
				,x).frobeniusInner(shapeFunction1.averageInDerivative(face,x))*(double)weight.value(x),
				face.cell1Ds,
				face.flatDimension,
				face.otherCoordinate);
		}
		if (name.equals(VALUE_NORMALAVERAGE_VALUE_NORMALAVERAGE))
		{
			return TPFaceIntegral.integrateNonTensorProduct(x -> shapeFunction1.normalAverageInValue(face
				,x).frobeniusInner(shapeFunction2.normalAverageInValue(face,x))*(double)weight.value(x),
				face.cell1Ds,
				face.flatDimension,
				face.otherCoordinate);
		}
		if (name.equals(BOUNDARY_VALUE))
		{
			if (face.isBoundaryFace())
				return TPFaceIntegral.integrateNonTensorProduct(x -> shapeFunction1.normalAverageInValue(face
					, x).frobeniusInner(shapeFunction2.normalAverageInValue(face, x)) * (double) weight.value(x),
					face.cell1Ds,
					face.flatDimension,
					face.otherCoordinate);
			else
				return 0;
		}
		throw new UnsupportedOperationException("unkown face integral name");
	}
	
}
