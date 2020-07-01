package linalg;

import java.util.stream.IntStream;

public class BlockMatrix
{
	protected DoubleTensor [][] blocks;
	protected final int blockN;
	protected final int blockSize;
	public BlockMatrix(int blockN, int blockSize)
	{
		this.blockSize =blockSize;
		this.blockN = blockN;
		blocks = new DoubleTensor[blockN][blockN];
	}
	
	public DoubleTensor[][] getBlocks()
	{
		return blocks;
	}
	
	public int getBlockN()
	{
		return blockN;
	}
	
	public int getBlockSize()
	{
		return blockSize;
	}
	
	public void initializeBlocks()
	{
		for(int i = 0; i < blockN; i++)
		{
			for(int j = 0; j < blockN; j++)
				blocks[i][j] = new DoubleTensor(blockSize,blockSize,true);
		}
	}
	public void copyFromSparse(DoubleTensor matrix)
	{

		IntStream.range(0,matrix.sparseEntries).parallel().forEach(i->{
			int blockY = matrix.sparseYs[i]/blockSize;
			int blockX = matrix.sparseXs[i]/blockSize;
			int subY = matrix.sparseYs[i]%blockSize;
			int subX = matrix.sparseXs[i]%blockSize;
			blocks[blockY][blockX].add(subY,subX,matrix.sparseValues[i]);
		});
	}
	public void copyFromDense(DoubleTensor matrix)
	{
		IntStream.range(0,matrix.getM()).parallel().forEach(i->{
			for(int j = 0; j < matrix.getN(); j++)
			{
				int blockY = i/blockSize;
				int blockX = j/blockSize;
				int subY = i%blockSize;
				int subX = j%blockSize;
				blocks[blockY][blockX].add(subY,subX,matrix.at(i,j));
			}
		});

	}
	public BlockMatrix(DoubleTensor matrix, int blockSize)
	{
		if(matrix.getDimension() != 2 || matrix.getM() != matrix.getN() || matrix.getN()%blockSize != 0)
			throw new UnsupportedOperationException();
		this.blockSize = blockSize;
		blockN = matrix.getN()/blockSize;
		blocks = new DoubleTensor[blockN][blockN];
		initializeBlocks();
		if(matrix.isSparse())
		{
			copyFromSparse(matrix);
		}
		else
		{
			copyFromDense(matrix);
		}

	}
}
