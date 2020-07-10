package tensorproduct;

import basic.*;
import linalg.CoordinateVector;

import java.util.List;
import java.util.function.ToDoubleFunction;

public class TPFaceIntegral extends FaceIntegral<TPCell<TPShapeFunction>,TPFace<TPShapeFunction>,TPShapeFunction>
{
	public static String VALUE_JUMP_VALUE_JUMP = "ValueJumpValueJump";
	public static String GRAD_NORMALAVERAGE_VALUE_JUMP = "GradNormalaverageValueJump";
	public static String VALUE_JUMP_GRAD_NORMALAVERAGE = "ValueJumpGradNormalaverage";
	public static String GRAD_VALUE_NORMAL = "GradValueNormal";
	public static String VALUE_GRAD_NORMAL = "ValueGradNormal";
	public static String VALUE_VALUE = "ValueValue";
	private final boolean weightIsTensorProduct;
	public TPFaceIntegral(Function<?,?,?> weight, String name, boolean weightIsTensorProduct)
	{
		super(weight,name);
		this.weightIsTensorProduct = weightIsTensorProduct;
		if(name.equals(VALUE_VALUE) && !(weight.value(new CoordinateVector(weight.getDomainDimension())) instanceof Double))
			throw new IllegalArgumentException();
		if(name.equals(VALUE_JUMP_GRAD_NORMALAVERAGE) && !(weight.value(new CoordinateVector(weight.getDomainDimension())) instanceof Double))
			throw new IllegalArgumentException();
		if(name.equals(GRAD_NORMALAVERAGE_VALUE_JUMP) && !(weight.value(new CoordinateVector(weight.getDomainDimension())) instanceof Double))
			throw new IllegalArgumentException();
		if(name.equals(GRAD_VALUE_NORMAL) && !(weight.value(new CoordinateVector(weight.getDomainDimension())) instanceof Double))
			throw new IllegalArgumentException();
		if(name.equals(VALUE_GRAD_NORMAL) && !(weight.value(new CoordinateVector(weight.getDomainDimension())) instanceof Double))
			throw new IllegalArgumentException();
		if(name.equals(VALUE_VALUE) && !(weight.value(new CoordinateVector(weight.getDomainDimension())) instanceof Double))
			throw new IllegalArgumentException();
	}
	public TPFaceIntegral(String name)
	{
		super(name);
		this.weightIsTensorProduct = true;
		if(name.equals(VALUE_VALUE) && !(weight.value(new CoordinateVector(weight.getDomainDimension())) instanceof Double))
			throw new IllegalArgumentException();
	}
	static double integrateTensorProduct(ToDoubleFunction<CoordinateVector> eval, List<Cell1D> cells,
	                                     int flatDimension, double otherCoordinate)
	{
		return TPCellIntegral.integrateTensorProduct(x->{
			double [] point = new double [cells.size()+1];
			int subd = 0;
			for (int j = 0; j < cells.size()+1; j++)
			{
				if(j == flatDimension)
					point[j] = otherCoordinate;
				else
					point[j] = x.at(subd++);
			}
			return eval.applyAsDouble(CoordinateVector.fromValues(point));
		},cells);
	}
	static double integrateNonTensorProduct(ToDoubleFunction<CoordinateVector> eval, List<Cell1D> cells,
	                                        int flatDimension, double otherCoordinate)
	{
		return TPCellIntegral.integrateNonTensorProduct(x->{
			double [] point = new double [cells.size()+1];
			int subd = 0;
			for (int j = 0; j < cells.size()+1; j++)
			{
				if(j == flatDimension)
					point[j] = otherCoordinate;
				else
					point[j] = x.at(subd++);
			}
			return eval.applyAsDouble(CoordinateVector.fromValues(point));
		},cells);
	}
	
	@Override
	public double evaluateFaceIntegral(TPFace<TPShapeFunction> face, TPShapeFunction shapeFunction1, TPShapeFunction shapeFunction2)
	{
		if(name.equals(VALUE_JUMP_VALUE_JUMP))
		{
			if(weightIsTensorProduct)
				return integrateTensorProduct(x->shapeFunction1.jumpInValue(face,x)*shapeFunction2.jumpInValue(face,x)*(Double) weight.value(x),
					face.cell1Ds,
					face.flatDimension,
					face.otherCoordinate);
			else
				return integrateNonTensorProduct(x->shapeFunction1.jumpInValue(face,x)*shapeFunction2.jumpInValue(face,x)*(Double) weight.value(x),
					face.cell1Ds,
					face.flatDimension,
					face.otherCoordinate);
		}
		if(name.equals(GRAD_NORMALAVERAGE_VALUE_JUMP))
		{
			if(weightIsTensorProduct)
				return integrateTensorProduct(x->shapeFunction1.normalAverageInDerivative(face,x)*shapeFunction2.jumpInValue(face,x)*(Double) weight.value(x),
					face.cell1Ds,
					face.flatDimension,
					face.otherCoordinate);
			else
				return integrateNonTensorProduct(x->shapeFunction1.normalAverageInDerivative(face,x)*shapeFunction2.jumpInValue(face,x)*(Double) weight.value(x),
					face.cell1Ds,
					face.flatDimension,
					face.otherCoordinate);
		}
		if(name.equals(VALUE_JUMP_GRAD_NORMALAVERAGE))
		{
			if(weightIsTensorProduct)
				return integrateTensorProduct(x->shapeFunction1.jumpInValue(face,x)*shapeFunction2.normalAverageInDerivative(face,x)*(Double) weight.value(x),
					face.cell1Ds,
					face.flatDimension,
					face.otherCoordinate);
			else
				return integrateNonTensorProduct(x->shapeFunction1.jumpInValue(face,x)*shapeFunction2.normalAverageInDerivative(face,x)*(Double) weight.value(x),
					face.cell1Ds,
					face.flatDimension,
					face.otherCoordinate);
		}
		if(name.equals(GRAD_VALUE_NORMAL))
		{
			if(weightIsTensorProduct)
				return integrateTensorProduct(x->shapeFunction1.gradient(x).inner(face.getNormal().value(x))*shapeFunction2.value(x)*(Double) weight.value(x),
					face.cell1Ds,
					face.flatDimension,
					face.otherCoordinate);
			else
				return integrateNonTensorProduct(x->shapeFunction1.gradient(x).inner(face.getNormal().value(x))*shapeFunction2.value(x)*(Double) weight.value(x),
					face.cell1Ds,
					face.flatDimension,
					face.otherCoordinate);
		}
		if(name.equals(VALUE_GRAD_NORMAL))
		{
			if (weightIsTensorProduct)
				return integrateTensorProduct(x -> shapeFunction2.gradient(x).inner(face.getNormal().value(x)) * shapeFunction1.value(x) * (Double) weight.value(x),
					face.cell1Ds,
					face.flatDimension,
					face.otherCoordinate);
			else
				return integrateNonTensorProduct(x -> shapeFunction2.gradient(x).inner(face.getNormal().value(x)) * shapeFunction1.value(x) * (Double) weight.value(x),
					face.cell1Ds,
					face.flatDimension,
					face.otherCoordinate);
		}
		if(name.equals(VALUE_VALUE))
		{
			if (weightIsTensorProduct)
				return integrateTensorProduct(x -> shapeFunction1.value(x) * shapeFunction1.value(x) * (Double) weight.value(x),
					face.cell1Ds,
					face.flatDimension,
					face.otherCoordinate);
			else
				return integrateNonTensorProduct(x ->
						shapeFunction1.value(x)* shapeFunction1.value(x) * (Double) weight.value(x),
					face.cell1Ds,
					face.flatDimension,
					face.otherCoordinate);
		}
		throw new UnsupportedOperationException("unkown face integral name");
	}

}
