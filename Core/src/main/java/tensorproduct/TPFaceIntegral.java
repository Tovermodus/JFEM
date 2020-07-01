package tensorproduct;

import basic.*;
import linalg.DoubleTensor;

public class TPFaceIntegral extends FaceIntegral
{
	TPFunction tpWeight;
	TPVectorFunction tpVectorWeight;
	public static String VALUE_JUMP_VALUE_JUMP = "ValueJumpValueJump";
	public static String GRAD_NORMALAVERAGE_VALUE_JUMP = "GradNormalaverageValueJump";
	public static String VALUE_JUMP_GRAD_NORMALAVERAGE = "ValueJumpGradNormalaverage";
	public static String GRAD_VALUE_NORMAL = "GradValueNormal";
	public static String VALUE_GRAD_NORMAL = "ValueGradNormal";
	public static String VALUE_VALUE = "ValueValue";
	public  TPFaceIntegral(ScalarFunction weight, String name)
	{
		super(weight, name);
	}
	public  TPFaceIntegral(TensorFunction vectorWeight, String name)
	{
		super(vectorWeight, name);
	}
	public TPFaceIntegral(TPFunction weight, String name)
	{
		super(weight, name);
		this.tpWeight = weight;
	}

	public TPFaceIntegral(TPVectorFunction vectorWeight, String name)
	{
		super(vectorWeight, name);
		this.tpVectorWeight = vectorWeight;
	}

	public TPFaceIntegral(String name)
	{
		this(new TPFunction(Function1D.oneFunction(),Function1D.oneFunction()), name);

	}

	@Override
	public double evaluateFaceIntegral(Face face, ScalarShapeFunction function1, ScalarShapeFunction function2)
	{
		if(face instanceof TPFace)
		{
			return this.evaluateFaceIntegral((TPFace)face, function1, function2);
		}
		throw new UnsupportedOperationException();
	}

	public double evaluateFaceIntegral(TPFace face, ScalarShapeFunction function1, ScalarShapeFunction function2)
	{
		Cell1D cell1d = face.getCell1d();
		int normaldirection = face.getNormaldirection();
		double otherCoordinate = face.getOtherCoordinate();
		if(name.equals(VALUE_JUMP_VALUE_JUMP))
		{

			double ret = 0;
			for(int i = 0; i < cell1d.weights.length; i++)
			{
				double weight = cell1d.weights[i];
				DoubleTensor point = new DoubleTensor(2);
				point.set(1-normaldirection,cell1d.points[i]);
				point.set(normaldirection,otherCoordinate);
				ret += face.jumpInValue(function1,point)*face.jumpInValue(function2,point)*this.weight.value(point)*weight;
			}
			return ret;
		}
		if(name.equals(GRAD_NORMALAVERAGE_VALUE_JUMP))
		{
			double ret = 0;
			for(int i = 0; i < cell1d.weights.length; i++)
			{
				double weight = cell1d.weights[i];
				DoubleTensor point = new DoubleTensor(2);
				point.set(1-normaldirection,cell1d.points[i]);
				point.set(normaldirection,otherCoordinate);
				ret += face.jumpInValue(function2,point)*face.normalAverageInDerivative(function1,
					point)*this.weight.value(point)*weight;
			}
			return ret;
		}
		if(name.equals(VALUE_JUMP_GRAD_NORMALAVERAGE))
		{

			double ret = 0;
			for(int i = 0; i < cell1d.weights.length; i++)
			{
				double weight = cell1d.weights[i];
				DoubleTensor point = new DoubleTensor(2);
				point.set(1-normaldirection,cell1d.points[i]);
				point.set(normaldirection,otherCoordinate);
				ret += face.jumpInValue(function1,point)*face.normalAverageInDerivative(function2,
					point)*this.weight.value(point)*weight;
			}
			return ret;
		}
		if(name.equals(GRAD_VALUE_NORMAL))
		{

			if(!face.isBoundaryFace())
				return 0;
			double ret = 0;
			for(int i = 0; i < cell1d.weights.length; i++)
			{
				double weight = cell1d.weights[i];
				DoubleTensor point = new DoubleTensor(2);
				point.set(1-normaldirection,cell1d.points[i]);
				point.set(normaldirection,otherCoordinate);
				ret += function1.derivative(point).inner(face.getNormal().value(point))*function2.value(point)*this.weight.value(point)*weight;
			}
			return ret;
		}
		if(name.equals(VALUE_GRAD_NORMAL))
		{

			if(!face.isBoundaryFace())
				return 0;
			double ret = 0;
			for(int i = 0; i < cell1d.weights.length; i++)
			{
				double weight = cell1d.weights[i];
				DoubleTensor point = new DoubleTensor(2);
				point.set(1-normaldirection,cell1d.points[i]);
				point.set(normaldirection,otherCoordinate);
				ret += function2.derivative(point).inner(face.getNormal().value(point))*function1.value(point)*this.weight.value(point)*weight;
			}
			return ret;
		}
		if(name.equals(VALUE_VALUE))
		{

			if(!face.isBoundaryFace())
				return 0;
			double ret = 0;
			for(int i = 0; i < cell1d.weights.length; i++)
			{
				double weight = cell1d.weights[i];
				DoubleTensor point = new DoubleTensor(2);
				point.set(1-normaldirection,cell1d.points[i]);
				point.set(normaldirection,otherCoordinate);
				ret += function2.value(point)*function1.value(point)*this.weight.value(point)*weight;
			}
			return ret;
		}
		throw new UnsupportedOperationException();
	}

	double evaluateFaceIntegral(Face face, TPShapeFunction function1,
	                            TPShapeFunction function2)
	{
		/*
		 *TODO!!!
		 */
		return evaluateFaceIntegral(face,(ScalarShapeFunction)function1,(ScalarShapeFunction)function2);
	}
}
