package tensorproduct;

import basic.*;
import linalg.CoordinateMatrix;
import linalg.CoordinateVector;
import linalg.Matrix;
import linalg.Vector;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;
import java.util.function.ToDoubleFunction;

public class TPCellIntegralViaReferenceCell<ST extends ScalarShapeFunctionWithReferenceShapeFunction<TPCell,TPFace,TPEdge,ST>> extends TPCellIntegral<ST>
{
	Map<ReferenceCellIdentificationTriplet<TPCell, TPFace, TPEdge, ST>, Double> savedValues;
	
	public TPCellIntegralViaReferenceCell(double weight, String name, boolean weightIsTensorProduct)
	{
		super(ScalarFunction.constantFunction(weight), name, weightIsTensorProduct);
		savedValues = new ConcurrentHashMap<>();
	}
	
	public TPCellIntegralViaReferenceCell(String name)
	{
		super(name);
		savedValues = new ConcurrentHashMap<>();
	}
	
	@Override
	public double evaluateCellIntegral(TPCell cell, ST shapeFunction1,
	                                   ST shapeFunction2)
	{
		ST referenceShapeFunction1 = shapeFunction1.getReferenceShapeFunctionRelativeTo(cell);
		ST referenceShapeFunction2 = shapeFunction2.getReferenceShapeFunctionRelativeTo(cell);
		ReferenceCellIdentificationTriplet<TPCell, TPFace, TPEdge, ST> key =
			new ReferenceCellIdentificationTriplet<>(referenceShapeFunction1,
				referenceShapeFunction2,
				cell.getReferenceCell());
		if (savedValues.containsKey(key))
			return savedValues.get(key);
		else
		{
			savedValues.put(key, super.evaluateCellIntegral(cell, shapeFunction1, shapeFunction2));
			return savedValues.get(key);
		}
	}
}