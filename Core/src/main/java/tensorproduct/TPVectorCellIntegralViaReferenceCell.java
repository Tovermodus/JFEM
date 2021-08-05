package tensorproduct;

import basic.*;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;

public class TPVectorCellIntegralViaReferenceCell<ST extends VectorShapeFunctionWithReferenceShapeFunction<TPCell, TPFace>> extends TPVectorCellIntegral<ST>
{
	Map<ReferenceCellIdentificationTriplet<TPCell, TPFace, VectorShapeFunctionWithReferenceShapeFunction<TPCell,TPFace>>,
		Double> savedValues;
	int i = 0;
	public TPVectorCellIntegralViaReferenceCell(double weight, String name)
	{
		super(ScalarFunction.constantFunction(weight), name);
		savedValues = new ConcurrentHashMap<>();
	}
	
	public TPVectorCellIntegralViaReferenceCell(String name)
	{
		super(name);
		savedValues = new ConcurrentHashMap<>();
	}
	
	@Override
	public double evaluateCellIntegral(TPCell cell, ST shapeFunction1,
	                                   ST shapeFunction2)
	{
		VectorShapeFunctionWithReferenceShapeFunction<TPCell,TPFace> referenceShapeFunction1
			= shapeFunction1.createReferenceShapeFunctionRelativeTo(cell);
		VectorShapeFunctionWithReferenceShapeFunction<TPCell,TPFace> referenceShapeFunction2
			= shapeFunction2.createReferenceShapeFunctionRelativeTo(cell);
		ReferenceCellIdentificationTriplet<TPCell, TPFace,
			VectorShapeFunctionWithReferenceShapeFunction<TPCell,TPFace>> key =
			new ReferenceCellIdentificationTriplet<>(referenceShapeFunction1,
				referenceShapeFunction2,
				cell.getReferenceCell());
		if (savedValues.containsKey(key))
		{
			//System.out.println(i++);
			return savedValues.get(key);
		}
		else
		{
			//System.out.println(i++);
			savedValues.put(key, super.evaluateCellIntegral(cell, shapeFunction1, shapeFunction2));
			return savedValues.get(key);
		}
	}
}
