package basic;

import linalg.CoordinateMatrix;
import linalg.CoordinateVector;
import linalg.Vector;

import java.util.HashMap;
import java.util.Map;

public class VectorFESpaceFunction<ST extends VectorShapeFunction<?, ?>> implements VectorFunction
{
	protected final HashMap<ST, Double> coefficients;
	
	public VectorFESpaceFunction(final ST[] functions, final double[] coefficients)
	{
		assert (functions.length == coefficients.length);
		this.coefficients = new HashMap<>();
		for (int i = 0; i < functions.length; i++)
		{
			this.coefficients.put(functions[i], coefficients[i]);
		}
	}
	
	public VectorFESpaceFunction(final Map<Integer, ST> functions, final Vector coefficients)
	{
		assert (functions.size() == coefficients.size());
		this.coefficients = new HashMap<>();
		for (final Map.Entry<Integer, ST> function : functions.entrySet())
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
	public CoordinateVector value(final CoordinateVector pos)
	{
		return coefficients
			.entrySet()
			.stream()
			.parallel()
			.map(entry -> entry.getKey().value(pos).mul(entry.getValue()))
			.reduce(new CoordinateVector(getDomainDimension()), CoordinateVector::add);
	}
	
	@Override
	public CoordinateMatrix gradient(final CoordinateVector pos)
	{
		return coefficients
			.entrySet()
			.stream()
			.parallel()
			.map(entry -> entry.getKey().gradient(pos).mul(entry.getValue()))
			.reduce(new CoordinateMatrix(getRangeDimension(), getDomainDimension()), CoordinateMatrix::add);
	}
	
	@Override
	public int getRangeDimension()
	{
		return coefficients.keySet().iterator().next().getRangeDimension();
	}
}
