package tensorproduct;

import basic.*;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;

public class TPFaceIntegralViaReferenceFace<ST extends ScalarShapeFunctionWithReferenceShapeFunction<TPCell, TPFace
	>> extends TPFaceIntegral<ST>
{
	Map<ReferenceFaceIdentificationTriplet<TPCell, TPFace, ScalarShapeFunctionWithReferenceShapeFunction<TPCell,TPFace
		>>, Double> savedValues;
	public TPFaceIntegralViaReferenceFace(double weight, String name)
	{
		super(ScalarFunction.constantFunction(weight),name);
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
		ScalarShapeFunctionWithReferenceShapeFunction<TPCell,TPFace
			> referenceShapeFunction1 = shapeFunction1.createReferenceShapeFunctionRelativeTo(face);
		ScalarShapeFunctionWithReferenceShapeFunction<TPCell,TPFace
			> referenceShapeFunction2 = shapeFunction2.createReferenceShapeFunctionRelativeTo(face);
		ReferenceFaceIdentificationTriplet<TPCell, TPFace, ScalarShapeFunctionWithReferenceShapeFunction<TPCell,TPFace
			>> key =
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
