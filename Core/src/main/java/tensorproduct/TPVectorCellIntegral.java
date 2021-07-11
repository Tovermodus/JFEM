package tensorproduct;

import basic.CellIntegral;
import basic.Function;
import basic.VectorShapeFunction;
import com.google.common.base.Stopwatch;
import linalg.CoordinateMatrix;
import linalg.CoordinateVector;

public class TPVectorCellIntegral<ST extends VectorShapeFunction<TPCell,TPFace,TPEdge,ST>> extends CellIntegral<TPCell,
	ST>
{
	public static final String GRAD_GRAD = "GradGrad";
	public static final String VALUE_VALUE = "ValueValue";
	public static final String ROT_ROT = "RotRot";
	public TPVectorCellIntegral(Function<?,?,?> weight, String name)
	{
		super(weight, name);
		if(name.equals(GRAD_GRAD) && !(weight.defaultValue() instanceof Double))
			throw new IllegalArgumentException();
		if(name.equals(VALUE_VALUE) && !(weight.defaultValue() instanceof Double))
			throw new IllegalArgumentException();
		if(name.equals(ROT_ROT) && !(weight.defaultValue() instanceof Double))
			throw new IllegalArgumentException();
	}
	public TPVectorCellIntegral(String name)
	{
		super(name);
	}
	@Override
	public double evaluateCellIntegral(TPCell cell, ST shapeFunction1, ST shapeFunction2)
	{
		if(name.equals(GRAD_GRAD))
		{
			return TPCellIntegral.integrateNonTensorProduct(x->shapeFunction1.gradient(x).frobeniusInner(shapeFunction2.gradient(x))*(Double)weight.value(x),cell.cell1Ds);
		}
		if(name.equals(VALUE_VALUE))
		{
			return TPCellIntegral.integrateNonTensorProduct(x->shapeFunction1.value(x).inner(shapeFunction2.value(x))*(Double)weight.value(x),cell.cell1Ds);
		}
		if(name.equals(ROT_ROT))
		{
			return TPCellIntegral.integrateNonTensorProduct(x->shapeFunction1.curl(x).inner(shapeFunction2.curl(x))*(Double)weight.value(x),cell.cell1Ds);
		}
		throw new UnsupportedOperationException("unknown integral name");
	}
}
