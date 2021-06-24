package tensorproduct;

import basic.*;
import linalg.CoordinateMatrix;
import linalg.CoordinateVector;
import linalg.Matrix;
import linalg.Vector;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.function.ToDoubleFunction;

public class TPCellIntegralViaReferenceCell<ST extends ScalarShapeFunctionWithReferenceShapeFunction<TPCell,TPFace,TPEdge,ST>> extends TPCellIntegral<ST>
{
	Map<ReferenceCellIdentificationTriplet<TPCell, TPFace, TPEdge, ST>, Double> savedValues;
	public TPCellIntegralViaReferenceCell(Function<?, ?, ?> weight, String name, boolean weightIsTensorProduct)
	{
		super(weight, name, weightIsTensorProduct);
		savedValues = new HashMap<>();
	}
	
	public TPCellIntegralViaReferenceCell(String name)
	{
		super(name);
		savedValues = new HashMap<>();
	}
	@Override
	public double evaluateCellIntegral(TPCell cell, ST shapeFunction1,
	                                   ST shapeFunction2)
	{
		if(name.equals(GRAD_GRAD))
		{
			ReferenceCellIdentificationTriplet<TPCell,TPFace,TPEdge, ST> key =
				new ReferenceCellIdentificationTriplet<>(shapeFunction1.getReferenceShapeFunctionRelativeTo(cell),
					shapeFunction2.getReferenceShapeFunctionRelativeTo(cell),
					cell.getReferenceCell());
			if(savedValues.containsKey(key))
				return savedValues.get(key);
			else{
				double conversionFactor = cell.getGradGradConversionFactor();
				TPCellIntegral<ST> refIntegral = new TPCellIntegral<ST>(weight, name, weightIsTensorProduct);
				savedValues.put(key, refIntegral.evaluateCellIntegral(cell, shapeFunction1,
					shapeFunction2)*conversionFactor);
			}
		}
		return super.evaluateCellIntegral(cell,shapeFunction1,shapeFunction2);
	}
}