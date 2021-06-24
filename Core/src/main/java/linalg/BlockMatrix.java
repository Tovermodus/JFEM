package linalg;

import basic.PerformanceArguments;
import com.google.common.primitives.Ints;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Stream;

public class BlockMatrix implements MutableMatrix
{
	DenseMatrix[][] blocks;
	
	public BlockMatrix(int blockSizeX, int blockSizeY, int blockNX, int blockNY)
	{
		blocks = new DenseMatrix[blockNX][blockNY];
		for (int i = 0; i < blockNX; i++)
			for (int j = 0; j < blockNY; j++)
				blocks[i][j] = new DenseMatrix(blockSizeX, blockSizeY);
	}
	
	public BlockMatrix(int blockSizeX, int blockSizeY, int blockNX, int blockNY, Matrix matrix)
	{
		this(blockSizeX, blockSizeY, blockNX, blockNY);
		copyFromMatrix(matrix);
	}
	
	public BlockMatrix(BlockMatrix matrix)
	{
		this(matrix.getBlockSizeX(), matrix.getBlockSizeY(), matrix.getBlockNX(), matrix.getBlockNY(), matrix);
	}
	
	public int getBlockSizeX()
	{
		return blocks[0][0].getRows();
	}
	
	public int getBlockSizeY()
	{
		return blocks[0][0].getCols();
	}
	
	public int getBlockNX()
	{
		return blocks.length;
	}
	
	public int getBlockNY()
	{
		return blocks[0].length;
	}
	
	private int[] blockCoords(int... coordinates)
	{
		int blockI = coordinates[0] / getBlockSizeX();
		int blockJ = coordinates[1] / getBlockSizeY();
		int subI = coordinates[0] % getBlockSizeX();
		int subJ = coordinates[1] % getBlockSizeY();
		return new int[]{blockI, blockJ, subI, subJ};
	}
	
	protected void copyFromMatrix(Matrix matrix)
	{
		for (int i = 0; i < matrix.getRows(); i++)
			for (int j = 0; j < matrix.getCols(); j++)
			{
				int[] cs = blockCoords(i, j);
				blocks[cs[0]][cs[1]].set(matrix.at(i, j), cs[2], cs[3]);
			}
	}
	
	@Override
	public double at(int... coordinates)
	{
		if (PerformanceArguments.getInstance().executeChecks)
			if (coordinates.length != 2)
				throw new IllegalArgumentException("Wrong number of coordinates");
		int[] cs = blockCoords(coordinates);
		return blocks[cs[0]][cs[1]].at(cs[2], cs[3]);
	}
	
	@Override
	public void set(double value, int... coordinates)
	{
		if (PerformanceArguments.getInstance().executeChecks)
			if (coordinates.length != 2)
				throw new IllegalArgumentException("Wrong number of coordinates");
		int[] cs = blockCoords(coordinates);
		blocks[cs[0]][cs[1]].set(value, cs[2], cs[3]);
	}
	
	@Override
	public void add(double value, int... coordinates)
	{
		if (PerformanceArguments.getInstance().executeChecks)
			if (coordinates.length != 2)
				throw new IllegalArgumentException("Wrong number of coordinates");
		int[] cs = blockCoords(coordinates);
		blocks[cs[0]][cs[1]].add(value, cs[2], cs[3]);
	}
	
	@Override
	public void addInPlace(Tensor other)
	{
		this.blocks = this.add(other).blocks;
	}
	
	@Override
	public void mulInPlace(double scalar)
	{
		Arrays.stream(blocks).parallel().forEach(matrixrow -> {
			for(DenseMatrix d: matrixrow)
				d.mulInPlace(scalar);
		});
	}
	
	@Override
	public BlockMatrix add(Tensor other)
	{
		if (PerformanceArguments.getInstance().executeChecks)
			if (!getShape().equals(other.getShape()))
				throw new IllegalArgumentException("Incompatible sizes");
		BlockMatrix ret = new BlockMatrix(this);
		for (int i = 0; i < ret.getRows(); i++)
			for (int j = 0; j < ret.getCols(); j++)
				ret.add(other.at(i, j), i, j);
		return ret;
	}
	
	@Override
	public BlockMatrix sub(Tensor other)
	{
		
		if (PerformanceArguments.getInstance().executeChecks)
			if (getShape() != other.getShape())
				throw new IllegalArgumentException("Incompatible sizes");
		BlockMatrix ret = new BlockMatrix(this);
		for (int i = 0; i < ret.getRows(); i++)
			for (int j = 0; j < ret.getCols(); j++)
				ret.add(-other.at(i, j), i, j);
		return ret;
	}
	
	@Override
	public BlockMatrix mul(double scalar)
	{
		BlockMatrix ret = new BlockMatrix(getBlockSizeX(), getBlockSizeY(), getBlockNX(), getBlockNY());
		for (int i = 0; i < ret.getRows(); i++)
			for (int j = 0; j < ret.getCols(); j++)
				ret.set(scalar * at(i, j), i, j);
		return ret;
	}
	
	@Override
	public List<Vector> unfoldDimension(int dimension)
	{
		if (PerformanceArguments.getInstance().executeChecks)
			if (dimension > 1)
				throw new IllegalArgumentException("Matrix is two dimensional");
		List<Vector> ret = new ArrayList<>();
		if (dimension == 0)
		{
			for (int i = 0; i < getRows(); i++)
			{
				DenseVector d = new DenseVector(getCols());
				for (int j = 0; j < getCols(); j++)
					d.add(at(i, j), j);
				ret.add(d);
			}
		} else
		{
			for (int i = 0; i < getCols(); i++)
			{
				DenseVector d = new DenseVector(getCols());
				for (int j = 0; j < getCols(); j++)
					d.add(at(j, i), j);
				ret.add(d);
			}
		}
		return ret;
	}
	
	@Override
	public int getSparseEntryCount()
	{
		return 0;
	}
	
	@Override
	public boolean isSparse()
	{
		return false;
	}
	
	@Override
	public IntCoordinates getShape()
	{
		return new IntCoordinates(getBlockNX() * getBlockSizeX(), getBlockNY() * getBlockSizeY());
	}
	
	@Override
	public long size()
	{
		return (long) getBlockSizeX() * getBlockSizeY() * getBlockNX() * getBlockNY();
	}
	
	
	@Override
	public Matrix transpose()
	{
		BlockMatrix ret = new BlockMatrix(getBlockSizeX(), getBlockSizeY(), getBlockNX(), getBlockNY());
		for (int i = 0; i < ret.getRows(); i++)
			for (int j = 0; j < ret.getCols(); j++)
				ret.set(at(j, i), i, j);
		return null;
	}
	
	@Override
	public Vector mvMul(Vector vector)
	{
		if (PerformanceArguments.getInstance().executeChecks)
			if (getShape().get(1) != vector.getShape().get(0))
				throw new IllegalArgumentException("Incompatible sizes");
		DenseVector ret = new DenseVector(getShape().get(0));
		for (int i = 0; i < getRows(); i++)
			for (int k = 0; k < getCols(); k++)
				ret.add(at(i, k) * vector.at(k), i);
		return ret;
	}
	
	@Override
	public Vector tvMul(Vector vector)
	{
		if (PerformanceArguments.getInstance().executeChecks) if (getShape().get(0) != vector.getShape().get(0))
			throw new IllegalArgumentException("Incompatible sizes");
		DenseVector ret = new DenseVector(getShape().get(1));
		for (int i = 0; i < getCols(); i++)
			for (int k = 0; k < getRows(); k++)
				ret.add(at(k, i) * vector.at(k), i);
		return ret;
	}
	
	@Override
	public Matrix mmMul(Matrix matrix)
	{
		if (PerformanceArguments.getInstance().executeChecks)
			if (getShape().get(1) != matrix.getShape().get(0))
				throw new IllegalArgumentException("Incompatible sizes");
		BlockMatrix ret = new BlockMatrix(getBlockSizeX(), getBlockSizeY(), getBlockNX(), getBlockNY());
		for (int i = 0; i < getRows(); i++)
			for (int j = 0; j < matrix.getCols(); j++)
				for (int k = 0; k < getCols(); k++)
					ret.add(at(i, k) * matrix.at(k, j), i, j);
		return ret;
	}
}
