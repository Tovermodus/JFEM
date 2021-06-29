package linalg;

public interface MutableVector extends Vector, MutableTensor
{
	default void setAll(double [] values, int startingPoint)
	{
		setAll(DenseVector.vectorFromValues(values), startingPoint);
	}
	
	default void setAll(double [] values, IntCoordinates startingPoint)
	{
		setAll(DenseVector.vectorFromValues(values), startingPoint.get(0));
	}
	
	default void setAll(Vector values, IntCoordinates startingPoint)
	{
		setAll(values, startingPoint.get(0));
	}
	default void setAll(Vector values, int startingPoint)
	{
		for (int i = 0; i < values.getLength(); i++)
		{
			set(values.at(i), startingPoint);
		}
	}
	default void addAll(double [] values, int startingPoint)
	{
		addAll(DenseVector.vectorFromValues(values), startingPoint);
	}
	
	default void addAll(double [] values, IntCoordinates startingPoint)
	{
		addAll(DenseVector.vectorFromValues(values), startingPoint.get(0));
	}
	
	default void addAll(Vector values, IntCoordinates startingPoint)
	{
		addAll(values, startingPoint.get(0));
	}
	default void addAll(Vector values, int startingPoint)
	{
		for (int i = 0; i < values.getLength(); i++)
		{
			add(values.at(i), startingPoint);
		}
	}
}
