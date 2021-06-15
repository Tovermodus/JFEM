package linalg;

import com.google.common.primitives.Ints;

import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.OptionalDouble;
import java.util.stream.IntStream;

public interface MutableTensor extends Tensor
{
	void set(double value, int... coordinates);
	void add(double value, int...coordinates);
	
	MutableTensor add(MutableTensor other);
	void addInPlace(MutableTensor other);
	void subInPlace(MutableTensor other);
	void mulInPlace(double scalar);
	default MutableTensor sub(MutableTensor other)
	{
		if(Ints.toArray(getShape()) != Ints.toArray(other.getShape()))
			throw new IllegalArgumentException("Tensors are of different size");
		return add(other.mul(-1.));
	}
	MutableTensor mul(double scalar);
	
}