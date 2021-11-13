package linalg;

public interface AddableVectorMultiplyable<M extends AddableVectorMultiplyable<M>>
	extends VectorMultiplyable
{
	M add(final M other);
	
	@Override
	M mul(final double scalar);
}
