package linalg;

public interface MutableVector
	extends Vector, MutableTensor
{
	default void setAll(final double[] values, final int startingPoint)
	{
		setAll(DenseVector.vectorFromValues(values), startingPoint);
	}
	
	default void setAll(final double[] values, final IntCoordinates startingPoint)
	{
		setAll(DenseVector.vectorFromValues(values), startingPoint.get(0));
	}
	
	default void setAll(final Vector values, final IntCoordinates startingPoint)
	{
		setAll(values, startingPoint.get(0));
	}
	
	default void setAll(final Vector values, final int startingPoint)
	{
		for (int i = 0; i < values.getLength(); i++)
		{
			set(values.at(i), startingPoint);
		}
	}
	
	default void addAll(final double[] values, final int startingPoint)
	{
		addAll(DenseVector.vectorFromValues(values), startingPoint);
	}
	
	default void addAll(final double[] values, final IntCoordinates startingPoint)
	{
		addAll(DenseVector.vectorFromValues(values), startingPoint.get(0));
	}
	
	default void addAll(final Vector values, final IntCoordinates startingPoint)
	{
		addAll(values, startingPoint.get(0));
	}
	
	default void addAll(final Vector values, final int startingPoint)
	{
		for (int i = 0; i < values.getLength(); i++)
		{
			add(values.at(i), startingPoint);
		}
	}
}
