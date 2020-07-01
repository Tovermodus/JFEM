package tensorproduct;
import basic.BoundaryFaceIntegral;
import basic.Face;
import basic.ScalarFunction;
import basic.ScalarShapeFunction;
import linalg.DoubleTensor;
import tensorproduct.TPFace;

public class TPBoundaryFaceIntegral extends BoundaryFaceIntegral
{
	public TPBoundaryFaceIntegral(ScalarFunction rightHandSide, ScalarFunction weight)
	{
		super(rightHandSide, weight);
	}

	@Override
	public double evaluateBoundaryFaceIntegral(Face f, ScalarShapeFunction v)
	{
		if(f instanceof TPFace)
		{
			System.out.println("make nice TBOUNFARYFACE");
			return this.evaluateBoundaryFaceIntegral((TPFace) f, v);
		}
		throw new UnsupportedOperationException();
	}

	public double evaluateBoundaryFaceIntegral(TPFace f, ScalarShapeFunction v)
	{
		if(!f.isBoundaryFace())
			return 0;
		double ret = 0;
		for(int i = 0; i < f.getCell1d().weights.length; i++)
		{
			double weight = f.getCell1d().weights[i];
			DoubleTensor point = new DoubleTensor(2);
			point.set(f.getNormaldirection(), f.getCell1d().points[i]);
			point.set(1-f.getNormaldirection(),f.getOtherCoordinate());
			ret += v.value(point)*rightHandSide.value(point)*this.weight.value(point)*weight;
		}
		return ret;
	}
}
