package linalg;

import java.util.stream.IntStream;

public class IsoBlockDiagonalMatrix extends IsoBlockMatrix
{
	public IsoBlockDiagonalMatrix(final int blockSize, final int blockN)
	{
		super(blockSize, blockSize, blockN, blockN);
	}
	
	public IsoBlockDiagonalMatrix(final int blockSize, final int blockN, final Matrix matrix)
	{
		this(blockSize, blockN);
		copyFromMatrix(matrix);
	}
	
	public IsoBlockDiagonalMatrix(final IsoBlockMatrix matrix)
	{
		this(matrix.getBlockSizeX(), matrix.getBlockNX(), matrix);
	}
	
	@Override
	public Vector mvMul(final Vector vector)
	{
		final BlockVector vect = new BlockVector(getBlockSizeX(), getBlockNX(), vector);
		return BlockVector.fromBlocks((DenseVector[]) IntStream
			.range(0, getBlockNX())
			.parallel()
			.mapToObj(i -> blocks[i][i].mvMul(vect.getBlocks()[i]))
			.toArray());
	}
	
	public Vector solve(final Vector vector)
	{
		final BlockVector vect = new BlockVector(getBlockSizeX(), getBlockNX(), vector);
		return BlockVector.fromBlocks((DenseVector[]) IntStream
			.range(0, getBlockNX())
			.parallel()
			.mapToObj(i -> blocks[i][i].solve(vect.getBlocks()[i]))
			.toArray());
	}
	
	public IsoBlockDiagonalMatrix inverse()
	{
		final IsoBlockDiagonalMatrix ret = new IsoBlockDiagonalMatrix(this);
		IntStream.range(0, getBlockNX()).parallel().forEach(i ->
		                                                    {
			                                                    ret.blocks[i][i] = blocks[i][i].inverse();
		                                                    });
		return ret;
	}
}
