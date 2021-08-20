package tensorproduct;

import basic.CellIntegral;
import basic.Function;
import basic.VectorShapeFunction;
import linalg.CoordinateMatrix;
import linalg.CoordinateVector;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

public class TPVectorCellIntegral<ST extends VectorShapeFunction<TPCell, TPFace>> extends CellIntegral<TPCell,
	ST>
{
	public static final String GRAD_GRAD = "GradGrad";
	public static final String VALUE_VALUE = "ValueValue";
	public static final String ROT_ROT = "RotRot";
	public static final String GRAD_VALUE = "GradValue";
	public static final String SYM_GRAD = "SymGrad";
	
	public TPVectorCellIntegral(final Function<?, ?, ?> weight, final String name)
	{
		this(weight, name, QuadratureRule1D.Gauss5);
	}
	
	public TPVectorCellIntegral(final String name)
	{
		this(name, QuadratureRule1D.Gauss5);
	}
	
	public TPVectorCellIntegral(final Function<?, ?, ?> weight, final String name, final QuadratureRule1D quadratureRule1D)
	{
		super(weight, name, quadratureRule1D);
		if (name.equals(GRAD_GRAD) && !(weight.defaultValue() instanceof Double))
			throw new IllegalArgumentException();
		if (name.equals(VALUE_VALUE) && !(weight.defaultValue() instanceof Double))
			throw new IllegalArgumentException();
		if (name.equals(ROT_ROT) && !(weight.defaultValue() instanceof Double))
			throw new IllegalArgumentException();
		if (name.equals(GRAD_VALUE) && !(weight.defaultValue() instanceof CoordinateVector))
			throw new IllegalArgumentException();
	}
	
	public TPVectorCellIntegral(final String name, final QuadratureRule1D quadratureRule1D)
	{
		super(name, quadratureRule1D);
	}
	
	@Override
	public double evaluateCellIntegral(final TPCell cell, final ST shapeFunction1, final ST shapeFunction2)
	{
		if (name.equals(SYM_GRAD))
		{
			return TPCellIntegral
				.integrateNonTensorProduct(x ->
				                           {
					                           final CoordinateMatrix grad1 =
						                           shapeFunction1.gradient(x);
					                           final CoordinateMatrix grad2 =
						                           shapeFunction2.gradient(x);
					
					                           return grad1
						                           .add(grad1.transpose())
						                           .frobeniusInner(grad2.add(grad2.transpose()))
						                           * (Double) weight.value(x) / 4;
				                           },
				                           cell,
				                           quadratureRule1D);
		}
		if (name.equals(GRAD_GRAD))
		{
			return TPCellIntegral.integrateNonTensorProduct(x -> shapeFunction1
				                                                .gradient(x)
				                                                .frobeniusInner(shapeFunction2.gradient(x))
				                                                * (Double) weight.value(x),
			                                                cell,
			                                                quadratureRule1D);
		}
		if (name.equals(VALUE_VALUE))
		{
			return TPCellIntegral.integrateNonTensorProduct(x -> shapeFunction1
				                                                .value(x)
				                                                .inner(shapeFunction2.value(x))
				                                                * (Double) weight.value(x),
			                                                cell,
			                                                quadratureRule1D);
		}
		if (name.equals(ROT_ROT))
		{
			return TPCellIntegral.integrateNonTensorProduct(x -> shapeFunction1
				                                                .curl(x)
				                                                .inner(shapeFunction2.curl(x))
				                                                * (Double) weight.value(x),
			                                                cell,
			                                                quadratureRule1D);
		}
		if (name.equals(GRAD_VALUE))
		{
			return TPCellIntegral.integrateNonTensorProduct(x ->
				                                                shapeFunction1
					                                                .gradient(x)
					                                                .mvMul(shapeFunction2.value(x))
					                                                .inner((CoordinateVector) weight
						                                                .value(x)),
			                                                cell,
			                                                quadratureRule1D);
		}
		throw new UnsupportedOperationException("unknown integral name");
	}
}
