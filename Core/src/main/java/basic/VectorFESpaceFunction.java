package basic;

import linalg.CoordinateVector;
import linalg.Vector;

import java.util.HashMap;
import java.util.Map;

public class VectorFESpaceFunction<ST extends VectorShapeFunction<?,?,?,ST>> extends VectorFunction
{
	private HashMap<ST, Double> coefficients;
	public VectorFESpaceFunction(ST[] functions, double[] coefficients)
	{
		assert(functions.length == coefficients.length);
		this.coefficients = new HashMap<>();
		for(int i = 0; i < functions.length; i++)
		{
			this.coefficients.put(functions[i], coefficients[i]);
		}
	}
	public VectorFESpaceFunction(Map<Integer, ST> functions, Vector coefficients)
	{
		assert(functions.size() == coefficients.size());
		this.coefficients = new HashMap<>();
		for(Map.Entry<Integer,ST> function:functions.entrySet())
		{
			this.coefficients.put(function.getValue(), coefficients.at(function.getKey()));
		}
	}
	
	@Override
	public int getDomainDimension()
	{
		return coefficients.keySet().iterator().next().getDomainDimension();
	}
	
	@Override
	public CoordinateVector value(CoordinateVector pos)
	{
		return coefficients.entrySet().stream().parallel().map(entry->entry.getKey().value(pos).mul(entry.getValue())).reduce(new CoordinateVector(getDomainDimension()), CoordinateVector::add);
	}
}
