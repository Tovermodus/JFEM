package linalg;

import com.google.common.primitives.Ints;

import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.OptionalDouble;
import java.util.stream.IntStream;

public interface MutableTensor
{
	void set(double value, IntCoordinates coordinates);
	void add(double value, IntCoordinates coordinates);
	
	void addInPlace(MutableTensor other);
	void subInPlace(MutableTensor other);
	void mulInPlace(double scalar);
	
}