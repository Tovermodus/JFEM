package mixed;

import basic.*;
import linalg.CoordinateVector;
import linalg.Vector;
import tensorproduct.TPCell;
import tensorproduct.TPFace;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;

public class MixedFESpaceFunction<CT extends Cell<CT,FT>, FT extends Face<CT,FT>, PF extends ScalarShapeFunction<CT,FT,
	PF>,
	VF extends VectorShapeFunction<CT,FT,VF>> extends MixedFunction
{
	private HashMap<MixedShapeFunction<CT,FT,PF,VF>, Double> coefficients;
	Map<MixedShapeFunction<CT, FT, PF, VF>, Double> pressureCoefficients;
	Map<MixedShapeFunction<CT, FT, PF, VF>, Double> velocityCoefficients;
	public MixedFESpaceFunction(MixedShapeFunction<CT,FT,PF,VF>[] functions, double[] coefficients)
	{
		super();
		assert(functions.length == coefficients.length);
		this.coefficients = new HashMap<>();
		for(int i = 0; i < functions.length; i++)
		{
			this.coefficients.put(functions[i], coefficients[i]);
		}
		initializeFunctionSets();
	}
	public MixedFESpaceFunction(Map<Integer, MixedShapeFunction<CT,FT,PF,VF>> functions, Vector coefficients)
	{
		assert(functions.size() == coefficients.size());
		this.coefficients = new HashMap<>();
		for(Map.Entry<Integer,MixedShapeFunction<CT,FT,PF,VF>> function:functions.entrySet())
		{
			this.coefficients.put(function.getValue(), coefficients.at(function.getKey()));
		}
		initializeFunctionSets();
	}
	private void initializeFunctionSets()
	{
		pressureCoefficients = new ConcurrentHashMap<>();
		velocityCoefficients = new ConcurrentHashMap<>();
		for (MixedShapeFunction<CT, FT, PF, VF> shapeFunction : coefficients.keySet())
		{
			if (shapeFunction.isPressure())
				pressureCoefficients.put(shapeFunction, coefficients.get(shapeFunction));
			else
				velocityCoefficients.put(shapeFunction, coefficients.get(shapeFunction));
		}
	}
	
	@Override
	public int getDomainDimension()
	{
		return coefficients.keySet().iterator().next().getDomainDimension();
	}
	
	@Override
	public MixedValue value(CoordinateVector pos)
	{
		PressureValue pf = new PressureValue(pressureCoefficients.entrySet().stream().parallel().mapToDouble(entry->entry.getKey().value(pos).getPressure()*entry.getValue()).sum());
		VelocityValue vf =
			velocityCoefficients.entrySet().stream().parallel().map(entry->(VelocityValue)(entry.getKey().value(pos)).mul(entry.getValue())).reduce(new VelocityValue(getDomainDimension()), (a,b)->(VelocityValue) (a.add(b)));
		return pf.add(vf);
	}
}
