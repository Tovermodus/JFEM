package multigrid;

import linalg.*;

public class BlockGaussseidelSmoother
	implements Smoother
{
	BlockDenseMatrix diagonalInverse;
	BlockSparseMatrix upper;
	BlockSparseMatrix lower;
	final int runs;
	final double omega = 1;
	
	public BlockGaussseidelSmoother(final BlockSparseMatrix blockMatrix, final int runs)
	{
		this.runs = runs;
		diagonalInverse = blockMatrix.getInvertedDiagonalMatrix();
		upper = new BlockSparseMatrix(blockMatrix,
		                              intCoordinates -> intCoordinates.get(0) < intCoordinates.get(1));
		lower = new BlockSparseMatrix(blockMatrix,
		                              intCoordinates -> intCoordinates.get(0) > intCoordinates.get(1));
	}
	
	@Override
	public Vector smooth(final VectorMultiplyable Operator,
	                     final Vector rhs,
	                     Vector iterate,
	                     final boolean verbose,
	                     final String prefix)
	{
		if (verbose)
			System.out.println(prefix + " " + Operator.mvMul(iterate)
			                                          .sub(rhs)
			                                          .euclidianNorm());
		for (int k = 0; k < runs; k++)
		{
			final var subRhs = BlockSparseMatrix.partitionVector(rhs.sub(upper.mvMul(iterate)
			                                                                  .mul(omega)),
			                                                     diagonalInverse.getBlockStarts(),
			                                                     diagonalInverse.getBlockEnds());
			iterate = new DenseVector(rhs.getLength());
			final DenseVector[] subRet = new DenseVector[subRhs.size()];
			for (int i = 0; i < subRhs.size(); i++)
			{
				final DenseVector thisRhs = new DenseVector(subRhs.get(i));
				for (int j = 0; j < i; j++)
				{
					if (lower.getBlockMatrix(i, j) != null)
						thisRhs.subInPlace(lower.getBlockMatrix(i, j)
						                        .mvMul(subRet[j])
						                        .mul(omega));
				}
				subRet[i] = diagonalInverse.getBlockMatrix(i, i)
				                           .mvMul(thisRhs);
				((DenseVector) iterate).addSmallVectorAt(subRet[i], lower.getBlockStarts()[i]);
			}
			if (verbose)
				System.out.println(prefix + " " + Operator.mvMul(iterate)
				                                          .sub(rhs)
				                                          .euclidianNorm());
		}
		return iterate;
	}
}
