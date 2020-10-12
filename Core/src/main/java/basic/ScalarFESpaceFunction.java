package basic;

import linalg.CoordinateVector;
import linalg.Vector;

import java.util.*;

public class ScalarFESpaceFunction<ST extends ScalarShapeFunction<?,?,?,ST>> extends ScalarFunction
{
	private HashMap<ST, Double> coefficients;
	public ScalarFESpaceFunction(ST[] functions, double[] coefficients)
	{
		assert(functions.length == coefficients.length);
		this.coefficients = new HashMap<>();
		for(int i = 0; i < functions.length; i++)
		{
			this.coefficients.put(functions[i], coefficients[i]);
		}
	}
	public ScalarFESpaceFunction(Map<Integer, ST> functions, Vector coefficients)
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
	public Double value(CoordinateVector pos)
	{
		return coefficients.entrySet().stream().parallel().mapToDouble(entry->(Double)(entry.getKey().value(pos))*entry.getValue()).sum();
	}
}
