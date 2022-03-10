package linalg;

import it.unimi.dsi.fastutil.ints.Int2ObjectArrayMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectMap;
import scala.Function2;

import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public abstract class AbstractSchurSolver<T extends VectorMultiplyable>
	implements VectorMultiplyable
{
	BlockSparseMatrix blockMatrix;
	final Int2ObjectMap<DenseMatrix> diagonalInverses;
	final int[] blockStarts;
	final int[] blockEnds;
	final int[] blockSizes;
	
	public AbstractSchurSolver(final BlockSparseMatrix blockMatrix)
	{
		this.blockMatrix = blockMatrix;
		diagonalInverses = new Int2ObjectArrayMap<>();
		blockStarts = blockMatrix.getBlockStarts();
		blockEnds = blockMatrix.getBlockEnds();
		blockSizes = blockMatrix.getBlockSizes();
		IntStream.range(1, blockStarts.length)
		         .parallel()
		         .forEach(i ->
		                  {
			                  System.out.println("inverting diagonal schur block " + i);
			                  final DenseMatrix inv = getDiagonalBlock(i - 1).inverse();
			                  synchronized (this)
			                  {
				                  diagonalInverses.put(i - 1, inv);
				                  if (diagonalInverses.get(i - 1) == null)
					                  throw new IllegalArgumentException(
						                  "matrix does not contain necessary diagonal block");
				                  if (getTopBlock(i - 1) == null)
					                  throw new IllegalArgumentException(
						                  "matrix does not contain necessary top block");
				                  if (getLeftBlock(i - 1) == null)
					                  throw new IllegalArgumentException(
						                  "matrix does not contain necessary left block");
			                  }
		                  });
		System.out.println("abstract schur built");
	}
	
	public int getBlockSize()
	{
		return diagonalInverses.size() + 1;
	}
	
	public SparseMatrix getDiagonalBlock(final int i)
	{
		return blockMatrix.getBlocks()
		                  .get(new IntCoordinates(blockStarts[i + 1], blockStarts[i + 1]));
	}
	
	public DenseMatrix getDiagonalInverse(final int i)
	{
		return diagonalInverses.get(i);
	}
	
	public SparseMatrix getTopBlock(final int i)
	{
		return blockMatrix.getBlocks()
		                  .get(new IntCoordinates(0, blockStarts[i + 1]));
	}
	
	public SparseMatrix getLeftBlock(final int i)
	{
		return blockMatrix.getBlocks()
		                  .get(new IntCoordinates(blockStarts[i + 1], 0));
	}
	
	public SparseMatrix getSchurBlock()
	{
		return blockMatrix.getBlockMatrix(0, 0);
	}
	
	public abstract T schurMvMul();
	
	protected abstract Function2<T, Vector, Vector> solveSchur();
	
	@Override
	public int getVectorSize()
	{
		return blockMatrix.getVectorSize();
	}
	
	@Override
	public int getTVectorSize()
	{
		return blockMatrix.getTVectorSize();
	}
	
	@Override
	public Vector mvMul(final Vector vector)
	{
		final List<DenseVector> subVectors = BlockSparseMatrix.partitionVector(vector,
		                                                                       blockStarts,
		                                                                       blockEnds);
		final DenseVector schurVector = subVectors.get(0);
		subVectors.remove(0);
		IntStream.range(0, diagonalInverses.size())
		         .parallel()
		         .mapToObj(i -> getTopBlock(i)
			         .mvMul(getDiagonalInverse(i)
				                .mvMul(subVectors.get(i))))
		         .forEach(schurVector::subInPlace);
		final List<DenseVector> solvedSubVectors =
			IntStream.range(0, diagonalInverses.size())
			         .parallel()
			         .mapToObj(i -> getDiagonalInverse(i)
				         .mvMul(subVectors.get(i)))
			         .collect(
				         Collectors.toList());
		final DenseVector solvedSchurVector = new DenseVector(solveSchur().apply(schurMvMul(), schurVector));
		final List<DenseVector> newSubVectors
			= IntStream.range(0, diagonalInverses.size())
			           .parallel()
			           .mapToObj(i -> solvedSubVectors.get(i)
			                                          .sub(getDiagonalInverse(i)
				                                               .mvMul(getLeftBlock(i)
					                                                      .mvMul(solvedSchurVector))))
			           .collect(Collectors.toList());
		final DenseVector ret = new DenseVector(getVectorSize());
		ret.addSmallVectorAt(solvedSchurVector, 0);
		IntStream.range(0, diagonalInverses.size())
		         .forEach(i -> ret.addSmallVectorAt(newSubVectors.get(i),
		                                            blockStarts[i + 1]));
		return ret;
	}
	
	@Override
	public Vector tvMul(final Vector vector)
	{
		throw new UnsupportedOperationException("not implemented yet");
	}
	
	public void resetOffDiagonals(final BlockSparseMatrix systemMatrix)
	{
		this.blockMatrix = systemMatrix;
	}
}
