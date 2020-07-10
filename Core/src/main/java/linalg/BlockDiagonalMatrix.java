package linalg;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

public class BlockDiagonalMatrix extends BlockMatrix
{
	public BlockDiagonalMatrix(int blockSize, int blockN) {
		super(blockSize, blockSize, blockN, blockN);
	}
	public BlockDiagonalMatrix(int blockSize, int blockN, Matrix matrix)
	{
		this(blockSize, blockN);
		copyFromMatrix(matrix);
	}
	public BlockDiagonalMatrix(BlockMatrix matrix)
	{
		this(matrix.getBlockSizeX(), matrix.getBlockNX(), matrix);
	}

	public Vector mvMul(Vector vector)
	{
		BlockVector vect = new BlockVector(getBlockSizeX(), getBlockNX(), vector);
		return BlockVector.fromBlocks((DenseVector[]) IntStream.range(0,getBlockNX()).parallel().mapToObj(i->blocks[i][i].mvMul(vect.getBlocks()[i])).toArray());
	}
	public Vector solve(Vector vector)
	{
		BlockVector vect = new BlockVector(getBlockSizeX(), getBlockNX(), vector);
		return BlockVector.fromBlocks((DenseVector[]) IntStream.range(0,getBlockNX()).parallel().mapToObj(i->blocks[i][i].solve(vect.getBlocks()[i])).toArray());
	}
	public BlockDiagonalMatrix inverse()
	{
		BlockDiagonalMatrix ret = new BlockDiagonalMatrix(this);
		IntStream.range(0,getBlockNX()).parallel().forEach(i->{ret.blocks[i][i] = blocks[i][i].inverse();});
		return ret;
	}
}
