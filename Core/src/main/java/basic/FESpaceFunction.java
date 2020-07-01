package basic;

import linalg.DoubleTensor;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

public class FESpaceFunction extends ScalarFunction
{
	private Map<ScalarShapeFunction, Double> coefficients;
	public FESpaceFunction(ScalarShapeFunction functions [], double coefficients [])
	{
		assert(functions.length == coefficients.length);
		this.coefficients = new HashMap<>();
		for(int i = 0; i < functions.length; i++)
		{
			this.coefficients.put(functions[i], coefficients[i]);
		}
	}
	public FESpaceFunction(ArrayList<ScalarShapeFunction> functions, DoubleTensor coefficients)
	{
		assert(functions.size() == coefficients.size());
		this.coefficients = new HashMap<>();
		for(int i = 0; i < functions.size(); i++)
		{
			this.coefficients.put(functions.get(i), coefficients.at(i));
		}
	}
	@Override
	public double value(DoubleTensor pos)
	{
		return coefficients.entrySet().stream().parallel().mapToDouble(entry->entry.getKey().value(pos)*entry.getValue()).sum();
	}
	@Override
	public DoubleTensor derivative(DoubleTensor pos)
	{
		throw new UnsupportedOperationException();
	}
}
