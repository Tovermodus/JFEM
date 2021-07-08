package mixed;

import basic.*;
import linalg.CoordinateMatrix;
import linalg.CoordinateVector;
import linalg.Vector;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.concurrent.ConcurrentHashMap;

public class MixedFESpaceFunction<CT extends Cell<CT,FT,ET>, FT extends Face<CT,FT,ET>, ET extends Edge<CT,FT,ET>,
	PF extends ScalarShapeFunction<CT,FT,ET,
	PF>,
	VF extends VectorShapeFunction<CT,FT,ET,VF>> extends MixedFunction
{
	private HashMap<MixedShapeFunction<CT,FT,ET,PF,VF>, Double> coefficients;
	Map<MixedShapeFunction<CT, FT,ET, PF, VF>, Double> pressureCoefficients;
	Map<MixedShapeFunction<CT, FT,ET, PF, VF>, Double> velocityCoefficients;
	public MixedFESpaceFunction(MixedShapeFunction<CT,FT,ET,PF,VF>[] functions, double[] coefficients)
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
	public MixedFESpaceFunction(Map<Integer, MixedShapeFunction<CT,FT,ET,PF,VF>> functions, Vector coefficients)
	{
		assert(functions.size() == coefficients.size());
		this.coefficients = new HashMap<>();
		for(Map.Entry<Integer,MixedShapeFunction<CT,FT,ET,PF,VF>> function:functions.entrySet())
		{
			this.coefficients.put(function.getValue(), coefficients.at(function.getKey()));
		}
		initializeFunctionSets();
	}
	private void initializeFunctionSets()
	{
		pressureCoefficients = new ConcurrentHashMap<>();
		velocityCoefficients = new ConcurrentHashMap<>();
		for (MixedShapeFunction<CT, FT,ET, PF, VF> shapeFunction : coefficients.keySet())
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
	public ScalarFunction getPressureFunction()
	{
		MixedFESpaceFunction<CT, FT,ET, PF, VF> me = this;
		return new ScalarFunction()
		{
			@Override
			public int getDomainDimension()
			{
				return me.getDomainDimension();
			}
			
			@Override
			public Double value(CoordinateVector pos)
			{
				return me.value(pos).getPressure();
			}
			
			@Override
			public CoordinateVector gradient(CoordinateVector pos)
			{
				return me.gradient(pos).getPressureGradient();
			}
		};
	}
	
	@Override
	public VectorFunction getVelocityFunction()
	{
		MixedFESpaceFunction<CT, FT, ET,PF, VF> me = this;
		return new VectorFunction()
		{
			@Override
			public int getDomainDimension()
			{
				return me.getDomainDimension();
			}
			
			@Override
			public CoordinateVector value(CoordinateVector pos)
			{
				return me.value(pos).getVelocity();
			}
			@Override
			public CoordinateMatrix gradient(CoordinateVector pos)
			{
				return me.gradient(pos).getVelocityGradient();
			}
		};
	}
	
	@Override
	public MixedValue value(CoordinateVector pos)
	{
		PressureValue pf = new PressureValue(pressureCoefficients.entrySet().stream().parallel().mapToDouble(entry->entry.getKey().value(pos).getPressure()*entry.getValue()).sum());
		VelocityValue vf =
			velocityCoefficients.entrySet().stream().parallel().map(entry->(VelocityValue)(entry.getKey().value(pos)).mul(entry.getValue())).reduce(new VelocityValue(getDomainDimension()), (a,b)->(VelocityValue) (a.add(b)));
		return pf.add(vf);
	}
	@Override
	public MixedGradient gradient(CoordinateVector pos)
	{
		Optional<CoordinateVector> pressureGradient =
			pressureCoefficients.entrySet().stream().parallel().map(entry->entry.getKey().gradient(pos).getPressureGradient().mul(entry.getValue())).reduce(CoordinateVector::add);
		PressureGradient pf = new PressureGradient(pos.mul(0));
		if(pressureGradient.isPresent())
			pf = new PressureGradient(pressureGradient.get());
		Optional<CoordinateMatrix> velocityGradient =
			velocityCoefficients.entrySet().stream().parallel().map(entry->entry.getKey().gradient(pos).getVelocityGradient().mul(entry.getValue())).reduce(CoordinateMatrix::add);
		VelocityGradient vf = new VelocityGradient(pos.outer(pos).mul(0));
		if(velocityGradient.isPresent())
			vf = new VelocityGradient(velocityGradient.get());
		return pf.add(vf);
	}
	
}
