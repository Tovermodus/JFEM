package linalg;

import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public abstract class ImplicitSchurSolver
	extends AbstractSchurSolver<VectorMultiplyable>
	implements VectorMultiplyable
{
	
	public ImplicitSchurSolver(final BlockSparseMatrix blockMatrix)
	{
		super(blockMatrix);
	}
	
	@Override
	protected VectorMultiplyable schurMvMul()
	{
		return new VectorMultiplyable()
		{
			@Override
			public int getVectorSize()
			{
				return blockSizes[0];
			}
			
			@Override
			public int getTVectorSize()
			{
				return blockSizes[0];
			}
			
			@Override
			public Vector mvMul(final Vector vector)
			{
				final List<DenseVector> subVectors =
					IntStream.range(0, diagonalInverses.size())
					         .parallel()
					         .mapToObj(i -> getTopBlock(i).mvMul(
						                                      getDiagonalInverse(i).mvMul(
							                                      getLeftBlock(i).mvMul(vector)))
					                                      .mul(-1))
					         .collect(Collectors.toList());
				final DenseVector v = getSchurBlock().mvMul(vector);
				for (int i = 0; i < diagonalInverses.size(); i++)
					v.addInPlace(subVectors.get(i));
				return v;
			}
			
			@Override
			public Vector tvMul(final Vector vector)
			{
				throw new UnsupportedOperationException("not implemented yet");
			}
		};
	}
}
