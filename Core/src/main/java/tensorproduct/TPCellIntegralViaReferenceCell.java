package tensorproduct;

import basic.*;

import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;

public class TPCellIntegralViaReferenceCell<ST extends ScalarShapeFunctionWithReferenceShapeFunction<TPCell,TPFace,TPEdge>> extends TPCellIntegral<ST>
{
	Map<ReferenceCellIdentificationTriplet<TPCell, TPFace, TPEdge, ScalarShapeFunctionWithReferenceShapeFunction<TPCell,TPFace,TPEdge>>,
		Double> savedValues;
	
	public TPCellIntegralViaReferenceCell(double weight, String name)
	{
		super(ScalarFunction.constantFunction(weight), name);
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
		ScalarShapeFunctionWithReferenceShapeFunction<TPCell,TPFace,TPEdge> referenceShapeFunction1
			= shapeFunction1.createReferenceShapeFunctionRelativeTo(cell);
		ScalarShapeFunctionWithReferenceShapeFunction<TPCell,TPFace,TPEdge> referenceShapeFunction2
			= shapeFunction2.createReferenceShapeFunctionRelativeTo(cell);
		ReferenceCellIdentificationTriplet<TPCell, TPFace, TPEdge,
			ScalarShapeFunctionWithReferenceShapeFunction<TPCell,TPFace,TPEdge>> key =
			new ReferenceCellIdentificationTriplet<>(referenceShapeFunction1,
				referenceShapeFunction2,
				cell.getReferenceCell());
		if (savedValues.containsKey(key))
		{
			return savedValues.get(key);
		}
		else
		{
			savedValues.put(key, super.evaluateCellIntegral(cell, shapeFunction1, shapeFunction2));
			return savedValues.get(key);
		}
	}
}