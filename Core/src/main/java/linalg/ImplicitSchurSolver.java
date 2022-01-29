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
				final SubVectorCalculator calc = new SubVectorCalculator(vector);
				final Thread t = new Thread(calc);
				t.start();
				final DenseVector v = getSchurBlock().mvMul(vector);
				try
				{
					t.join();
				} catch (final InterruptedException e)
				{
					e.printStackTrace();
				}
				for (int i = 0; i < diagonalInverses.size(); i++)
					v.addInPlace(calc.getSubVectors()
					                 .get(i));
				return v;
			}
			
			@Override
			public Vector tvMul(final Vector vector)
			{
				throw new UnsupportedOperationException("not implemented yet");
			}
		};
	}
	
	class SubVectorCalculator
		implements Runnable
	{
		public List<DenseVector> subVectors;
		public Vector vector;
		
		public SubVectorCalculator(final Vector v)
		{
			vector = v;
		}
		
		public List<DenseVector> getSubVectors()
		{
			return subVectors;
		}
		
		@Override
		public void run()
		{
			subVectors = IntStream.range(0, diagonalInverses.size())
			                      .parallel()
			                      .mapToObj(i -> getTopBlock(i).mvMul(
				                                                   getDiagonalInverse(i).mvMul(
					                                                   getLeftBlock(i).mvMul(vector)))
			                                                   .mul(-1))
			                      .collect(Collectors.toList());
		}
	}
}
