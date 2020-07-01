package linalg;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

public class BlockDiagonalMatrix extends BlockMatrix implements VectorMultiplication
{
	public BlockDiagonalMatrix(int blockN, int blockSize) {super(blockN,blockSize);}
	public BlockDiagonalMatrix(DoubleTensor matrix, int blockSize)
	{
		super(matrix, blockSize);
	}

	public DoubleTensor mvmul(DoubleTensor vector)
	{
		if(vector.getDimension() != 1)
			throw new UnsupportedOperationException();
		ArrayList<DoubleTensor> blockVectors = new ArrayList<>();
		for(int i = 0; i < blockN; i++)
			blockVectors.add(new DoubleTensor(blockSize));
		IntStream.range(0,(int)vector.size()).forEach(i->
		{
			blockVectors.get((int) (i / blockSize)).add(i%blockSize,vector.at(i));
		});
		Stream<DoubleTensor> dtStream =
			IntStream.range(0,blockN).parallel().mapToObj(i->blocks[i][i].mvmul(blockVectors.get(i)));
		List<DoubleTensor> resultVectors = dtStream.collect(Collectors.toList());
		DoubleTensor retVector = new DoubleTensor(vector.size());
		IntStream.range(0,(int)vector.size()).parallel().forEach(i->
		{
			retVector.add(i, resultVectors.get((int) (i / blockSize)).at(i%blockSize));
		});
		return retVector;
	}
	public DoubleTensor solve(DoubleTensor vector)
	{
		if(vector.getDimension() != 1)
			throw new UnsupportedOperationException();
		ArrayList<DoubleTensor> blockVectors = new ArrayList<>();
		for(int i = 0; i < blockN; i++)
			blockVectors.add(new DoubleTensor(blockSize));
		IntStream.range(0,(int)vector.size()).forEach(i->
		{
			blockVectors.get((int) (i / blockSize)).add(i%blockSize,vector.at(i));
		});
		Stream<DoubleTensor> dtStream =
			IntStream.range(0,blockN).parallel().mapToObj(i->blocks[i][i].solve(blockVectors.get(i)));
		List<DoubleTensor> resultVectors = dtStream.collect(Collectors.toList());
		DoubleTensor retVector = new DoubleTensor(vector.size());
		IntStream.range(0,(int)vector.size()).parallel().forEach(i->
		{
			retVector.add(i, resultVectors.get((int) (i / blockSize)).at(i%blockSize));
		});
		return retVector;
	}
	public BlockDiagonalMatrix inverse()
	{
		BlockDiagonalMatrix ret = new BlockDiagonalMatrix(this.blockN,this.blockSize);
		System.out.println(this.blockSize);
		IntStream.range(0,blockN).parallel().forEach(i->{ret.blocks[i][i] = blocks[i][i].inverse();});
		return ret;
	}
	@Override
	public void initializeBlocks()
	{
		for(int i = 0; i < blockN; i++)
		{
			blocks[i][i] = new DoubleTensor(blockSize,blockSize,true);
		}
	}
	@Override
	public void copyFromSparse(DoubleTensor matrix)
	{

		IntStream.range(0,matrix.sparseEntries).parallel().forEach(i->{
			int blockY = matrix.sparseYs[i]/blockSize;
			int blockX = matrix.sparseXs[i]/blockSize;
			int subY = matrix.sparseYs[i]%blockSize;
			int subX = matrix.sparseXs[i]%blockSize;
			if(blockX == blockY)
				blocks[blockY][blockX].add(subY,subX,matrix.sparseValues[i]);
		});
	}
}
