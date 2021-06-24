package linalg;

import basic.PerformanceArguments;
import com.google.common.primitives.Ints;

import java.util.List;

public class BlockVector implements MutableVector
{
	DenseVector[] blocks;
	public BlockVector(int blockSize, int blockN)
	{
		blocks = new DenseVector[blockN];
		for (int i = 0; i < blockN; i++)
		{
			blocks[i] = new DenseVector(blockSize);
		}
	}
	public int getBlockSize()
	{
		return blocks[0].getLength();
	}
	public int getBlockN()
	{
		return blocks.length;
	}
	private int[] blockCoords(int...coordinates)
	{
		int blockI = coordinates[0] / getBlockSize();
		int subI = coordinates[0] % getBlockSize();
		return new int[]{blockI,subI};
	}
	public BlockVector(int blockSize, int blockN, Vector vect)
	{
		this(blockSize,blockN);
		copyFromVector(vect);
	}
	public BlockVector(BlockVector vect)
	{
		this(vect.getBlockSize(),vect.getBlockN(), vect);
		copyFromVector(vect);
	}
	private void copyFromVector(Vector vector)
	{
		for(int j = 0; j < vector.getLength(); j++)
		{
			int[] cs = blockCoords(j);
			blocks[cs[0]].set(vector.at(j),cs[1]);
		}
	}
	
	DenseVector[] getBlocks()
	{
		return blocks;
	}
	static BlockVector fromBlocks(DenseVector[] blocks)
	{
		if (PerformanceArguments.getInstance().executeChecks)
		for(int i = 1; i < blocks.length; i++)
			if(blocks[i].size() != blocks[0].size())
				throw new IllegalArgumentException("Blocks are not the same size");
		BlockVector ret = new BlockVector(blocks[0].getLength(),blocks.length);
		ret.blocks = blocks.clone();
		return ret;
	}
	
	@Override
	public int getSparseEntryCount()
	{
		return -1;
	}
	
	@Override
	public boolean isSparse()
	{
		return false;
	}
	
	@Override
	public IntCoordinates getShape()
	{
		return new IntCoordinates(getBlockN()*getBlockSize());
	}
	
	@Override
	public double at(int... coordinates)
	{
		if (PerformanceArguments.getInstance().executeChecks)
		if(coordinates.length != 1)
			throw new IllegalArgumentException("Wrong number of coordinates");
		int[] cs = blockCoords(coordinates);
		return blocks[cs[0]].at(cs[1]);
	}
	
	@Override
	public void set(double value, int... coordinates)
	{
		if (PerformanceArguments.getInstance().executeChecks)
		if(coordinates.length != 1)
			throw new IllegalArgumentException("Wrong number of coordinates");
		int[] cs = blockCoords(coordinates);
		blocks[cs[0]].set(value, cs[1]);
	}
	
	@Override
	public void add(double value, int... coordinates)
	{
		if (PerformanceArguments.getInstance().executeChecks)
		if(coordinates.length != 1)
			throw new IllegalArgumentException("Wrong number of coordinates");
		int[] cs = blockCoords(coordinates);
		blocks[cs[0]].add(value, cs[1]);
	}
	
	@Override
	public void addInPlace(Tensor other)
	{
		this.blocks = this.add(other).blocks;
	}
	
	@Override
	public void mulInPlace(double scalar)
	{
		this.blocks = this.mul(scalar).blocks;
	}
	
	@Override
	public BlockVector add(Tensor other)
	{
		if (PerformanceArguments.getInstance().executeChecks)
		if(getShape() != other.getShape())
			throw new IllegalArgumentException("Vectors are of different size");
		BlockVector ret = new BlockVector(this);
		if(!other.isSparse())
			for (int i = 0; i < getLength(); i++)
			{
				ret.add(other.at(i),i);
			}
		else
			for(IntCoordinates key: other.getCoordinateEntryList().keySet())
				ret.add(other.getCoordinateEntryList().get(key), key);
		return ret;
	}
	@Override
	public BlockVector sub(Tensor other)
	{
		if (PerformanceArguments.getInstance().executeChecks)
		if(getShape() != other.getShape())
			throw new IllegalArgumentException("Vectors are of different size");
		BlockVector ret = new BlockVector(this);
		if(!other.isSparse())
			for (int i = 0; i < getLength(); i++)
			{
				ret.add(-other.at(i),i);
			}
		else
			for(IntCoordinates key: other.getCoordinateEntryList().keySet())
				ret.add(-other.getCoordinateEntryList().get(key), key);
		return ret;
	}
	@Override
	public BlockVector mul(double scalar)
	{
		BlockVector ret = new BlockVector(getBlockSize(),getBlockN());
		for (int i = 0; i < getLength(); i++)
		{
			ret.set(at(i)*scalar,i);
		}
		return ret;
	}
	
	@Override
	public Matrix outer(Vector other)
	{
		throw new UnsupportedOperationException("outer not yet implemented");
	}
	
}
