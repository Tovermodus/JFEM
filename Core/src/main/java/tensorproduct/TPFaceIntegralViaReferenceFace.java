package tensorproduct;

import basic.*;

import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;

public class TPFaceIntegralViaReferenceFace<ST extends ScalarShapeFunctionWithReferenceShapeFunction<TPCell,TPFace,
	TPEdge>> extends TPFaceIntegral<ST>
{
	Map<ReferenceFaceIdentificationTriplet<TPCell, TPFace, TPEdge, ScalarShapeFunctionWithReferenceShapeFunction<TPCell,TPFace,
		TPEdge>>, Double> savedValues;
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
		ScalarShapeFunctionWithReferenceShapeFunction<TPCell,TPFace,
			TPEdge> referenceShapeFunction1 = shapeFunction1.createReferenceShapeFunctionRelativeTo(face);
		ScalarShapeFunctionWithReferenceShapeFunction<TPCell,TPFace,
			TPEdge> referenceShapeFunction2 = shapeFunction2.createReferenceShapeFunctionRelativeTo(face);
		ReferenceFaceIdentificationTriplet<TPCell, TPFace, TPEdge, ScalarShapeFunctionWithReferenceShapeFunction<TPCell,TPFace,
			TPEdge>> key =
			new ReferenceFaceIdentificationTriplet<>(referenceShapeFunction1,
				referenceShapeFunction2,
				face.getReferenceFace());
		if (savedValues.containsKey(key))
		{
			if(Math.abs(savedValues.get(key) -  super.evaluateFaceIntegral(face, shapeFunction1,
			shapeFunction2)) > 1e-5)
			{
				System.out.println();
				System.out.println();
				System.out.println();
				System.out.println(shapeFunction1.getGlobalIndex() + " " + shapeFunction2.getGlobalIndex());
				System.out.println(referenceShapeFunction1);
				System.out.println(shapeFunction1);
				System.out.println(referenceShapeFunction2);
				System.out.println(shapeFunction2);
				System.out.println(face);
				System.out.println(face.getReferenceFace());
				System.out.println(savedValues.get(key) + " " + super.evaluateFaceIntegral(face,
					shapeFunction1, shapeFunction2) + " integrals");
			}
			return savedValues.get(key);
		}
		else
		{
			savedValues.put(key, super.evaluateFaceIntegral(face, shapeFunction1, shapeFunction2));
			return savedValues.get(key);
		}
	}

}
