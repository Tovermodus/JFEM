package tensorproduct;

import basic.*;
import linalg.CoordinateVector;

import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;
import java.util.function.ToDoubleFunction;

public class TPFaceIntegralViaReferenceFace<ST extends ScalarShapeFunctionWithReferenceShapeFunction<TPCell,TPFace,
	TPEdge,ST>> extends TPFaceIntegral<ST>
{
	Map<ReferenceFaceIdentificationTriplet<TPCell, TPFace, TPEdge, ST>, Double> savedValues;
	public TPFaceIntegralViaReferenceFace(double weight, String name, boolean weightIsTensorProduct)
	{
		super(ScalarFunction.constantFunction(weight),name, weightIsTensorProduct);
		savedValues = new ConcurrentHashMap<>();
	}
	public TPFaceIntegralViaReferenceFace(String name)
	{
		super(name);
		savedValues = new ConcurrentHashMap<>();
	}
	
	@Override
	public double evaluateFaceIntegral(TPFace face, ST shapeFunction1,
	                                   ST shapeFunction2)
	{
		ST referenceShapeFunction1 = shapeFunction1.getReferenceShapeFunctionRelativeTo(face);
		ST referenceShapeFunction2 = shapeFunction2.getReferenceShapeFunctionRelativeTo(face);
		ReferenceFaceIdentificationTriplet<TPCell, TPFace, TPEdge, ST> key =
			new ReferenceFaceIdentificationTriplet<>(referenceShapeFunction1,
				referenceShapeFunction2,
				face.getReferenceFace());
		if (savedValues.containsKey(key))
		{
			return savedValues.get(key);
		}
		else
		{
			System.out.println(savedValues.size());
			savedValues.put(key, super.evaluateFaceIntegral(face, shapeFunction1, shapeFunction2));
			return savedValues.get(key);
		}
	}

}
