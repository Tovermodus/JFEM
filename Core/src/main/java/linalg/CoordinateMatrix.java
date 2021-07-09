package linalg;

import basic.PerformanceArguments;

import java.util.Arrays;

public class CoordinateMatrix extends DenseMatrix
{
	public CoordinateMatrix(int i, int j)
	{
		super(i, j);
		if (i > 3 || j > 3)
			throw new IllegalArgumentException("only 1D, 2D and 3D supported");
		
	}
	
	public CoordinateMatrix(Matrix matrix)
	{
		super(matrix);
		if (matrix.getShape().get(0) > 3 || matrix.getShape().get(1) > 3)
			throw new IllegalArgumentException("only 1D, 2D and 3D supported");
	}
	
	public CoordinateMatrix(double[][] matrix)
	{
		super(matrix);
	}
	
	public static CoordinateMatrix fromValues(int rows, int cols, double... vals)
	{
		if (PerformanceArguments.getInstance().executeChecks)
		if(rows*cols != vals.length)
			throw new IllegalArgumentException("matrix does not fit size");
		CoordinateMatrix ret = new CoordinateMatrix(rows, cols);
		for (int i = 0; i < rows * cols; i++)
		{
			ret.entries[i / cols][i % cols] = vals[i];
		}
		return ret;
	}
	
	@Override
	public CoordinateVector mvMul(Vector vector)
	{
		return (CoordinateVector) super.mvMul(vector);
	}
	
	@Override
	public CoordinateMatrix add(Tensor other)
	{
		if (!getShape().equals(other.getShape()))
			throw new IllegalArgumentException("Incompatible sizes");
		CoordinateMatrix ret = new CoordinateMatrix(this);
		for (int i = 0; i < ret.getRows(); i++)
			for (int j = 0; j < ret.getCols(); j++)
				ret.add(other.at(i, j), i, j);
		return ret;
	}
	
	@Override
	public CoordinateMatrix sub(Tensor other)
	{
		
		if (!getShape().equals(other.getShape()))
			throw new IllegalArgumentException("Incompatible sizes");
		CoordinateMatrix ret = new CoordinateMatrix(this);
		for (int i = 0; i < ret.getRows(); i++)
			for (int j = 0; j < ret.getCols(); j++)
				ret.add(-other.at(i, j), i, j);
		return ret;
	}
	
	public CoordinateMatrix inverse()
	{
		return new CoordinateMatrix(super.inverse());
	}
	
	@Override
	public CoordinateMatrix mul(double scalar)
	{
		CoordinateMatrix ret = new CoordinateMatrix(entries.length, entries[0].length);
		for (int i = 0; i < ret.getRows(); i++)
			for (int j = 0; j < ret.getCols(); j++)
				ret.set(scalar * at(i, j), i, j);
		return ret;
	}
	
	public double determinant()
	{
		if (this.getCols() == 1)
			return at(0, 0);
		if (this.getCols() == 2)
			return at(0, 0) * at(1, 1) - at(1, 0) * at(0, 1);
		if (this.getCols() == 3)
			return at(0, 0) * at(1, 1) * at(2, 2)
				+ at(0, 1) * at(1, 2) * at(2, 0)
				+ at(0, 2) * at(1, 0) * at(2, 1)
				- at(0, 0) * at(1, 2) * at(2, 1)
				- at(1, 0) * at(2, 2) * at(0, 1)
				- at(2, 0) * at(0, 2) * at(1, 1);
		throw new IllegalStateException("Dimension not allowed");
	}
}
