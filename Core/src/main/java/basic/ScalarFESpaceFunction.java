package basic;

import linalg.CoordinateVector;
import linalg.Vector;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class ScalarFESpaceFunction<ST extends ScalarShapeFunction<?,?,ST>> extends ScalarFunction
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
	public ScalarFESpaceFunction(List<ST> functions, Vector coefficients)
	{
		assert(functions.size() == coefficients.size());
		this.coefficients = new HashMap<>();
		for(int i = 0; i < functions.size(); i++)
		{
			this.coefficients.put(functions.get(i), coefficients.at(i));
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
