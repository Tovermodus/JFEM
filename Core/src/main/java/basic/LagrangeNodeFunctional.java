package basic;

import linalg.CoordinateMatrix;
import linalg.CoordinateVector;

import java.util.Objects;

public class LagrangeNodeFunctional implements NodeFunctional<Double, CoordinateVector, CoordinateMatrix>
{
	private final CoordinateVector point;
	
	public LagrangeNodeFunctional(final CoordinateVector point)
	{
		this.point = point;
	}
	
	public CoordinateVector getPoint()
	{
		return point;
	}
	
	@Override
	public FunctionSignature getFunctionSignature()
	{
		return new FunctionSignature(Double.class, CoordinateVector.class, CoordinateMatrix.class);
	}
	
	@Override
	public double evaluate(final Function<Double, CoordinateVector, CoordinateMatrix> func)
	{
		return func.value(point);
	}
	
	@Override
	public boolean usesFace(final Face<?, ?> f)
	{
		return f.isOnFace(point);
	}
	
	@Override
	public boolean equals(final Object o)
	{
		if (this == o) return true;
		if (o == null || getClass() != o.getClass()) return false;
		final LagrangeNodeFunctional that = (LagrangeNodeFunctional) o;
		return Objects.equals(point, that.point);
	}
	
	@Override
	public int hashCode()
	{
		return Objects.hash(point);
	}
	
	@Override
	public String toString()
	{
		return "LagrangeNodeFunctional{" +
			"point=" + point +
			'}';
	}
}
