package linalg;

import com.google.common.primitives.Ints;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

public class DenseMatrix implements Matrix, Decomposable, DirectlySolvable
{
	double[][] entries;
	private org.ujmp.core.Matrix ujmpmat = null;
	public DenseMatrix(int i, int j)
	{
		entries = new double[i][j];
	}
	public DenseMatrix(Matrix matrix)
	{
		entries = new double[matrix.getRows()][matrix.getCols()];
		if(matrix.isSparse())
		{
			for(Map.Entry<List<Integer>,Double> entry:matrix.getCoordinateEntryList().entrySet())
				entries[entry.getKey().get(0)][entry.getKey().get(1)] += entry.getValue();
		}
		else
		{
			for (int i = 0; i < getRows(); i++)
				for (int j = 0; j < getCols(); j++)
					entries[i][j] = matrix.at(i, j);
		}
	}
	public static DenseMatrix squareMatrixFromValues(double... values)
	{
		if(Math.pow((int)Math.sqrt(values.length),2) != values.length)
			throw new UnsupportedOperationException("nope");
		DenseMatrix ret = new DenseMatrix((int)Math.sqrt(values.length),(int)Math.sqrt(values.length));
		for(int i = 0; i < ret.getRows(); i++)
			for(int j = 0; j < ret.getCols(); j++)
				ret.set(values[i*ret.getRows()+j],i,j);
		return ret;
	}
	public static DenseMatrix identity(int n)
	{
		DenseMatrix ret = new DenseMatrix(n,n);
		for(int  i = 0; i < n; i++)
			ret.add(1,i,i);
		return ret;
	}
	@Override
	public double at(int... coordinates)
	{
		if(coordinates.length != 2)
			throw new IllegalArgumentException("Wrong number of coordinates");
		return entries[coordinates[0]][coordinates[1]];
	}
	
	@Override
	public synchronized void set(double value, int... coordinates)
	{
		if(coordinates.length != 2)
			throw new IllegalArgumentException("Wrong number of coordinates");
		entries[coordinates[0]][coordinates[1]] = value;
	}
	
	@Override
	public synchronized void add(double value, int... coordinates)
	{
		if(coordinates.length != 2)
			throw new IllegalArgumentException("Wrong number of coordinates");
		entries[coordinates[0]][coordinates[1]] += value;
	}
	
	@Override
	public DenseMatrix add(Tensor other)
	{
		if(!getShape().equals(other.getShape()))
			throw new IllegalArgumentException("Incompatible sizes");
		DenseMatrix ret = new DenseMatrix(this);
		for(int i = 0; i < ret.getRows(); i++)
			for(int j = 0; j < ret.getCols(); j++)
				ret.add(other.at(i,j),i,j);
		return ret;
	}
	
	@Override
	public DenseMatrix sub(Tensor other){
		
		if(!getShape().equals(other.getShape()))
			throw new IllegalArgumentException("Incompatible sizes");
		DenseMatrix ret = new DenseMatrix(this);
		for(int i = 0; i < ret.getRows(); i++)
			for(int j = 0; j < ret.getCols(); j++)
				ret.add(-other.at(i,j),i,j);
		return ret;
	}
	
	@Override
	public DenseMatrix mul(double scalar)
	{
		
		DenseMatrix ret = new DenseMatrix(entries.length, entries[0].length);
		for(int i = 0; i < ret.getRows(); i++)
			for(int j = 0; j < ret.getCols(); j++)
				ret.set(scalar*at(i,j),i,j);
		return ret;
	}
	
	@Override
	public List<Vector> unfoldDimension(int dimension)
	{
		if(dimension > 1)
			throw new IllegalArgumentException("Matrix is two dimensional");
		List<Vector> ret = new ArrayList<>();
		if(dimension == 0)
		{
			for (int i = 0; i < getRows(); i++)
			{
				ret.add(new DenseVector(getShape().get(1 )));
				for (int j = 0; j < getCols(); j++)
					ret.get(i).add(at(i,j),j);
			}
		}
		else
		{
			for (int i = 0; i < getCols(); i++)
			{
				ret.add(new DenseVector(getRows()));
				for (int j = 0; j < getRows(); j++)
					ret.get(i).add(at(j, i), j);
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
	public List<Integer> getShape()
	{
		return Ints.asList(entries.length, entries[0].length);
	}
	
	@Override
	public long size()
	{
		return entries.length*entries[0].length;
	}
	
	
	@Override
	public DenseMatrix transpose()
	{
		DenseMatrix ret = new DenseMatrix(entries[0].length, entries.length);
		for(int i = 0; i < ret.getRows(); i++)
			for(int j = 0; j < ret.getCols(); j++)
				ret.set(at(j,i),i,j);
		return ret;
	}
	
	@Override
	public DenseVector mvMul(Vector vector)
	{
		if(getCols()!=(vector.getLength()))
			throw new IllegalArgumentException("Incompatible sizes");
		DenseVector ret = new DenseVector(getRows());
		for(int i = 0; i < getRows(); i++)
			for(int k = 0; k < getCols(); k++)
				ret.add(at(i,k)*vector.at(k),i);
		return ret;
	}
	
	@Override
	public DenseVector tvMul(Vector vector)
	{
		if(getRows()!=(vector.getLength()))
			throw new IllegalArgumentException("Incompatible sizes");
		DenseVector ret = new DenseVector(getCols());
		for(int i = 0; i < getCols(); i++)
			for(int k = 0; k < getRows(); k++)
				ret.add(at(k,i)*vector.at(k),i);
		return ret;
	}
	
	@Override
	public DenseMatrix mmMul(Matrix matrix)
	{
		if(getCols()!=(matrix.getRows()))
			throw new IllegalArgumentException("Incompatible sizes");
		DenseMatrix ret = new DenseMatrix(getRows(), matrix.getCols());
		for(int i = 0; i < getRows(); i++)
			for(int j = 0; j < matrix.getCols(); j++)
				for(int k = 0; k < getCols(); k++)
					ret.add(at(i,k)*matrix.at(k,j),i,j);
		return ret;
	}
	
	
	@Override
	public DenseMatrix getLowerTriangleMatrix()
	{
		if(getRows() != getCols())
			throw new UnsupportedOperationException("cant get triangle from rectangle");
		DenseMatrix ret = new DenseMatrix(getRows(), getCols());
		for (int i = 0; i < getRows(); i++)
			for(int j = 0; j < i; j++)
				ret.add(at(i,j),i,j);
		return ret;
	}
	@Override
	public DenseMatrix getUpperTriangleMatrix()
	{
		if(getRows() != getCols())
			throw new UnsupportedOperationException("cant get triangle from rectangle");
		DenseMatrix ret = new DenseMatrix(getRows(), getCols());
		for (int i = 0; i < getRows(); i++)
			for(int j = i+1; j < getCols(); j++)
				ret.add(at(i,j),i,j);
		return ret;
	}
	
	@Override
	public DenseMatrix getDiagonalMatrix()
	{
		if(getRows() != getCols())
			throw new UnsupportedOperationException("cant get triangle from rectangle");
		DenseMatrix ret = new DenseMatrix(getRows(), getCols());
		for (int i = 0; i < getRows(); i++)
			ret.add(at(i,i),i,i);
		return ret;
	}
	private org.ujmp.core.Matrix toUJMPMatrix()
	{
		if(ujmpmat == null)
		{
			ujmpmat = org.ujmp.core.DenseMatrix.Factory.zeros(getRows(), getCols());
			for (int i = 0; i < getRows(); i++)
				for (int j = 0; j < getCols(); j++)
					ujmpmat.setAsDouble(at(i, j), i, j);
		}
		return ujmpmat;
	}
	private static DenseMatrix fromUJMPMatrix(org.ujmp.core.Matrix matrix)
	{
		DenseMatrix ret = new DenseMatrix((int)matrix.getRowCount(), (int)matrix.getColumnCount());
		
		for(int i = 0; i < ret.getRows(); i++)
			for(int j = 0; j < ret.getCols(); j++)
				ret.set(matrix.getAsDouble(i,j), i, j);
		return ret;
	}
	@Override
	public DenseMatrix inverse()
	{
		return fromUJMPMatrix(toUJMPMatrix().inv());
	}
	
	@Override
	public DenseVector solve(Vector rhs)
	{
		return DenseVector.fromUJMPVector(toUJMPMatrix().solve(new DenseVector(rhs).toUJMPVector()));
	}
	
	@Override
	public Vector solveSymm(Vector rhs)
	{
		return DenseVector.fromUJMPVector(toUJMPMatrix().solveSymm(new DenseVector(rhs).toUJMPVector()));
	}
	
	@Override
	public boolean equals(Object obj)
	{
		if(!(obj instanceof Matrix))
			return false;
		return almostEqual((Tensor) obj);
	}
	
	@Override
	public String toString()
	{
		return printFormatted();
	}
}
