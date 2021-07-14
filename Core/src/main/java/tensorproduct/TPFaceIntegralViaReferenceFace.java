package tensorproduct;

import basic.*;
import linalg.CoordinateVector;

import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;
import java.util.function.ToDoubleFunction;

public class TPFaceIntegralViaReferenceFace<ST extends ScalarShapeFunctionWithReferenceShapeFunction<TPCell,TPFace,
	TPEdge>> extends TPFaceIntegral<ST>
{
	Map<ReferenceFaceIdentificationTriplet<TPCell, TPFace, TPEdge, ScalarShapeFunctionWithReferenceShapeFunction<TPCell,TPFace,
		TPEdge>>, Double> savedValues;
	public TPFaceIntegralViaReferenceFace(double weight, String name, boolean weightIsTensorProduct)
	{
		super(ScalarFunction.constantFunction(weight),name, weightIsTensorProduct);
		savedValues = new ConcurrentHashMap<>();
	}
	public TPFaceIntegralViaReferenceFace(double weight, String name)
	{
		super(ScalarFunction.constantFunction(weight),name, true);
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
		ScalarShapeFunctionWithReferenceShapeFunction<TPCell,TPFace,
			TPEdge> referenceShapeFunction1 = shapeFunction1.getReferenceShapeFunctionRelativeTo(face);
		ScalarShapeFunctionWithReferenceShapeFunction<TPCell,TPFace,
			TPEdge> referenceShapeFunction2 = shapeFunction2.getReferenceShapeFunctionRelativeTo(face);
		ReferenceFaceIdentificationTriplet<TPCell, TPFace, TPEdge, ScalarShapeFunctionWithReferenceShapeFunction<TPCell,TPFace,
			TPEdge>> key =
			new ReferenceFaceIdentificationTriplet<>(referenceShapeFunction1,
				referenceShapeFunction2,
				face.getReferenceFace());
		if (savedValues.containsKey(key))
		{
			return savedValues.get(key);
		}
		else
		{
			savedValues.put(key, super.evaluateFaceIntegral(face, shapeFunction1, shapeFunction2));
			return savedValues.get(key);
		}
	}

}
