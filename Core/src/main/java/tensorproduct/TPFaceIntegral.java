package tensorproduct;

import basic.*;
import linalg.CoordinateVector;
import tensorproduct.geometry.Cell1D;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

import java.util.List;
import java.util.function.ToDoubleFunction;

public class TPFaceIntegral<ST extends ScalarShapeFunction<TPCell, TPFace>> extends FaceIntegral<TPFace,ST>
{
	public static String VALUE_JUMP_VALUE_JUMP = "ValueJumpValueJump";
	public static String INTERIOR_VALUE_JUMP_VALUE_JUMP = "InteriorValueJumpValueJump";
	public static String GRAD_NORMALAVERAGE_VALUE_JUMP = "GradNormalaverageValueJump";
	public static String VALUE_JUMP_GRAD_NORMALAVERAGE = "ValueJumpGradNormalaverage";
	public static String GRAD_VALUE_NORMAL = "GradValueNormal";
	public static String VALUE_GRAD_NORMAL = "ValueGradNormal";
	public static String VALUE_VALUE = "ValueValue";
	public static String BOUNDARY_VALUE = "BoundaryValue";
	
	public TPFaceIntegral(double weight, String name, QuadratureRule1D quadratureRule1D)
	{
		this(ScalarFunction.constantFunction(weight), name,quadratureRule1D);
	}
	public TPFaceIntegral(double weight, String name)
	{
		this(ScalarFunction.constantFunction(weight), name);
	}
	
	public TPFaceIntegral(Function<?, ?, ?> weight, String name)
	{
		this(weight, name, QuadratureRule1D.Gauss5);
	}
	
	public TPFaceIntegral(Function<?, ?, ?> weight, String name, QuadratureRule1D quadratureRule1D)
	{
		super(weight,name, quadratureRule1D);
		if(name.equals(VALUE_VALUE) && !(weight.value(new CoordinateVector(weight.getDomainDimension())) instanceof Double))
			throw new IllegalArgumentException();
		if(name.equals(BOUNDARY_VALUE) && !(weight.value(new CoordinateVector(weight.getDomainDimension())) instanceof Double))
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
		this(name, QuadratureRule1D.Gauss5);
	}
	
	public TPFaceIntegral(String name, QuadratureRule1D quadratureRule1D)
	{
		super(name, quadratureRule1D);
		if(name.equals(VALUE_VALUE) && !(weight.value(new CoordinateVector(weight.getDomainDimension())) instanceof Double))
			throw new IllegalArgumentException();
		if(name.equals(BOUNDARY_VALUE) && !(weight.value(new CoordinateVector(weight.getDomainDimension())) instanceof Double))
			throw new IllegalArgumentException();
	}
	public static double integrateNonTensorProduct(ToDoubleFunction<CoordinateVector> eval, TPFace f, QuadratureRule1D quadratureRule)
	{
		return integrateNonTensorProduct(eval, f.getComponentCells(), f.flatDimension, f.otherCoordinate,
			quadratureRule);
	}
	public static double integrateNonTensorProduct(ToDoubleFunction<CoordinateVector> eval, List<Cell1D> cells,
	                                               int flatDimension, double otherCoordinate,
	                                               QuadratureRule1D quadratureRule)
	{
		return TPCellIntegral.integrateNonTensorProduct(x ->
		{
			double[] point = new double[cells.size() + 1];
			int subd = 0;
			for (int j = 0; j < cells.size() + 1; j++)
			{
				if (j == flatDimension)
					point[j] = otherCoordinate;
				else
					point[j] = x.at(subd++);
			}
			return eval.applyAsDouble(CoordinateVector.fromValues(point));
		}, cells, quadratureRule);
	}
	
	@Override
	public double evaluateFaceIntegral(TPFace face, ST shapeFunction1,
	                                   ST shapeFunction2)
	{
		if (name.equals(VALUE_JUMP_VALUE_JUMP))
		{
			return integrateNonTensorProduct(x -> shapeFunction1.jumpInValue(face, x) * shapeFunction2.jumpInValue(face, x) * (Double) weight.value(x),
				face,
				quadratureRule1D);
		}
		if (name.equals(INTERIOR_VALUE_JUMP_VALUE_JUMP))
		{
			if(face.isBoundaryFace())
				return 0;
			return integrateNonTensorProduct(x -> shapeFunction1.jumpInValue(face, x) * shapeFunction2.jumpInValue(face, x) * (Double) weight.value(x),
				face,
				quadratureRule1D);
		}
		if (name.equals(GRAD_NORMALAVERAGE_VALUE_JUMP))
		{
			return integrateNonTensorProduct(x -> shapeFunction1.normalAverageInDerivative(face, x) * shapeFunction2.jumpInValue(face, x) * (Double) weight.value(x),
				face,
				quadratureRule1D);
		}
		if (name.equals(VALUE_JUMP_GRAD_NORMALAVERAGE))
		{
			return integrateNonTensorProduct(x -> shapeFunction1.jumpInValue(face, x) * shapeFunction2.normalAverageInDerivative(face, x) * (Double) weight.value(x),
				face,
				quadratureRule1D);
		}
		if (name.equals(GRAD_VALUE_NORMAL))
		{
			return integrateNonTensorProduct(x -> shapeFunction1.gradient(x).inner(face.getNormal().value(x)) * shapeFunction2.value(x) * (Double) weight.value(x),
				face,
				quadratureRule1D);
		}
		if (name.equals(VALUE_GRAD_NORMAL))
		{
			return integrateNonTensorProduct(x -> shapeFunction2.gradient(x).inner(face.getNormal().value(x)) * shapeFunction1.value(x) * (Double) weight.value(x),
				face,
				quadratureRule1D);
		}
		if (name.equals(VALUE_VALUE))
		{
			return integrateNonTensorProduct(x ->
					shapeFunction1.value(x) * shapeFunction2.value(x) * (Double) weight.value(x),
				face,
				quadratureRule1D);
		}
		if (name.equals(BOUNDARY_VALUE))
		{
			if(face.isBoundaryFace())
				return integrateNonTensorProduct(x ->
						shapeFunction1.jumpInValue(face, x)* shapeFunction2.jumpInValue(face,
							x)* (Double) weight.value(x),
					face,
					quadratureRule1D);
			return 0;
		}
		throw new UnsupportedOperationException("unkown face integral name");
	}

}
