package tensorproduct;

import basic.RightHandSideIntegral;
import basic.ScalarFunction;
import basic.ScalarShapeFunction;
import basic.Cell;
import com.google.common.collect.Table;
import linalg.DoubleTensor;

public class TPRightHandSideIntegral extends RightHandSideIntegral
{
	public TPRightHandSideIntegral(ScalarFunction rightHandSide, ScalarFunction weight)
	{
		super(rightHandSide, weight);
	}

	public TPRightHandSideIntegral(ScalarFunction rightHandSide)
	{
		super(rightHandSide);
	}

	@Override
	public double evaluateRightHandSideIntegral(Cell k, ScalarShapeFunction v)
	{
		if(k instanceof  TPCell)
			return this.evaluateRightHandSideIntegral((TPCell) k, v);
		throw new UnsupportedOperationException();
	}

	public double evaluateRightHandSideIntegral(TPCell k, ScalarShapeFunction v)
	{

		double ret = 0;
		DoubleTensor point;
		for(int i = 0; i < k.cellx.weights.length; i++)
			for(int j = 0; j < k.celly.weights.length; j++)
			{
				point = DoubleTensor.vectorFromValues(k.cellx.points[i],k.celly.points[j]);
				ret += v.value(point)*rightHandSide.value(point)*weight.value(point)*k.cellx.weights[i]*k.celly.weights[j];
			}
		return ret;
	}
}
