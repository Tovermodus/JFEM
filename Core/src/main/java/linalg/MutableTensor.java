package linalg;

import com.google.common.primitives.Ints;

import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.OptionalDouble;
import java.util.stream.IntStream;

public interface MutableTensor
{
	void set(double value, int... coordinates);
	void add(double value, int... coordinates);
	default void set(double value, IntCoordinates coordinates)
	{
		set(value, coordinates.asArray());
	}
	default void add(double value, IntCoordinates coordinates)
	{
		add(value, coordinates.asArray());
	}
	
	void addInPlace(Tensor other);
	default void subInPlace(Tensor other)
	{
		addInPlace(other.mul(-1.));
	}
	void mulInPlace(double scalar);
	
}