package linalg;

import com.google.common.primitives.Ints;

import java.util.*;

public class SparseMatrix implements Matrix, DirectlySolvable, Decomposable
{
	private final int rows;
	private final int cols;
	volatile double [] sparseValues;
	volatile int [] sparseXs;
	volatile int [] sparseYs;
	volatile int sparseEntries;
	public SparseMatrix(int rows, int cols)
	{
		this.rows = rows;
		this.cols = cols;
		sparseValues = new double[1];
		sparseYs = new int[1];
		sparseXs = new int[1];
		this.sparseEntries = 1;
	}
	public SparseMatrix(SparseMatrix m)
	{
		this.rows = m.rows;
		this.cols = m.cols;
		this.sparseEntries = m.sparseEntries;
		sparseXs = new int[sparseEntries];
		sparseYs = new int[sparseEntries];
		sparseValues = new double[sparseEntries];
		for(int i = 0;i < sparseEntries; i++)
		{
			sparseXs[i] = m.sparseXs[i];
			sparseYs[i] = m.sparseYs[i];
			sparseValues[i] = m.sparseValues[i];
		}
	}
	
	public static SparseMatrix identity(int i)
	{
		SparseMatrix ret = new SparseMatrix(i,i);
		for (int j = 0; j < i; j++)
		{
			ret.add(1,j,j);
		}
		return ret;
	}
	
	@Override
	public DenseMatrix inverse()
	{
		return new DenseMatrix(this).inverse();
	}
	
	@Override
	public Vector solve(Vector rhs)
	{
		return new DenseMatrix(this).solve(rhs);
	}
	
	@Override
	public Vector solveSymm(Vector rhs)
	{
		return new DenseMatrix(this).solveSymm(rhs);
	}
	
	@Override
	public double at(int... coordinates)
	{
		if(coordinates.length != 2)
			throw new IllegalArgumentException("Wrong number of coordinates");
		double ret = 0;
		for(int i = 0; i < sparseEntries; i++)
		{
			if(sparseYs[i] == coordinates[0])
				if(sparseXs[i] == coordinates[1])
					ret += sparseValues[i];
		}
		return ret;
	}
	
	private synchronized void resizeSparse()
	{
		double [] sparseVals = new double[sparseValues.length*2];
		int [] sparseX  = new int[sparseValues.length*2];
		int [] sparseY =  new int[sparseValues.length*2];
		for(int i = 0; i < sparseEntries; i++)
		{
			sparseVals[i] = sparseValues[i];
			sparseX[i] = sparseXs[i];
			sparseY[i] = sparseYs[i];
		}
		sparseValues = sparseVals;
		sparseYs = sparseY;
		sparseXs = sparseX;
	}
	
	@Override
	public synchronized void set(double value, int... coordinates)
	{
		if(coordinates.length != 2)
			throw new IllegalArgumentException("Wrong number of coordinates");
		if(sparseEntries >= sparseValues.length-2)
			resizeSparse();
		sparseValues[sparseEntries] = value - at(coordinates[0],coordinates[1]);
		sparseYs[sparseEntries] = coordinates[0];
		sparseXs[sparseEntries] = coordinates[1];
		sparseEntries++;
	}
	
	@Override
	public synchronized void add(double value, int... coordinates)
	{
		if(coordinates.length != 2)
			throw new IllegalArgumentException("Wrong number of coordinates");
		if(sparseEntries >= sparseValues.length-2)
			resizeSparse();
		sparseValues[sparseEntries] = value;
		sparseYs[sparseEntries] = coordinates[0];
		sparseXs[sparseEntries] = coordinates[1];
		sparseEntries++;
	}
	
	@Override
	public Matrix add(Tensor other)
	{
		if(!getShape().equals(other.getShape()))
			throw new IllegalArgumentException("Incompatible sizes");
		if(!(other instanceof SparseMatrix))
			return ((Matrix)other).add(this);
		SparseMatrix ret = new SparseMatrix(this);
		for (int i = 0; i < ((SparseMatrix) other).sparseEntries; i++)
		{
			ret.add(((SparseMatrix) other).sparseValues[i],((SparseMatrix) other).sparseYs[i],
				((SparseMatrix) other).sparseXs[i]);
		}
		return ret;
	}
	
	@Override
	public Map<List<Integer>, Double> getCoordinateEntryList()
	{
		Map<List<Integer>,Double> list = new HashMap<>();
		for (int i = 0; i < sparseEntries; i++)
		{
			list.put(Ints.asList(sparseYs[i],sparseXs[i]),sparseValues[i]);
		}
		return list;
	}
	
	@Override
	public SparseMatrix mul(double scalar)
	{
		SparseMatrix ret = new SparseMatrix(this);
		for (int i = 0; i < sparseEntries; i++)
		{
			ret.sparseValues[i] *= scalar;
		}
		return ret;
	}
	
	@Override
	public List<Vector> unfoldDimension(int dimension)
	{
		if(dimension > 1)
			throw new IllegalArgumentException("Matrix is two dimensional");
		List<Vector> ret = new ArrayList<>(getShape().get(dimension));
		for(int i = 0; i < getShape().get(dimension); i++)
			ret.add(new DenseVector(getShape().get(1-dimension)));
		if(dimension == 0)
			for (int i = 0; i < sparseEntries; i++)
				ret.get(sparseYs[i]).add(sparseValues[i],sparseXs[i]);
		if(dimension == 1)
			for (int i = 0; i < sparseEntries; i++)
				ret.get(sparseXs[i]).add(sparseValues[i],sparseYs[i]);
			
		return ret;
	}
	
	@Override
	public int getSparseEntryCount()
	{
		return sparseEntries;
	}
	
	@Override
	public boolean isSparse()
	{
		return true;
	}
	
	@Override
	public List<Integer> getShape()
	{
		return Ints.asList(rows,cols);
	}
	
	@Override
	public SparseMatrix transpose()
	{
		SparseMatrix ret = new SparseMatrix(cols, rows);
		ret.sparseEntries = sparseEntries;
		ret.sparseValues = sparseValues.clone();
		ret.sparseXs = sparseYs.clone();
		ret.sparseYs = sparseXs.clone();
		return ret;
	}
	
	@Override
	public Vector mvMul(Vector vector)
	{
		if(getCols()!=(vector.getLength()))
			throw new IllegalArgumentException("Incompatible sizes");
		DenseVector ret = new DenseVector(getRows());
		for (int i = 0; i < sparseEntries; i++)
			ret.add(sparseValues[i]*vector.at(sparseXs[i]),sparseYs[i]);
		return ret;
		
	}
	@Override
	public Vector tvMul(Vector vector)
	{
		if(getCols()!=(vector.getLength()))
			throw new IllegalArgumentException("Incompatible sizes");
		DenseVector ret = new DenseVector(getCols());
		for (int i = 0; i < sparseEntries; i++)
			ret.add(sparseValues[i]*vector.at(sparseYs[i]),sparseXs[i]);
		return ret;
		
	}
	
	@Override
	public Matrix mmMul(Matrix matrix)
	{
		if(getCols()!=(matrix.getRows()))
			throw new IllegalArgumentException("Incompatible sizes");
		if(!(matrix instanceof SparseMatrix))
			return matrix.tmMul(transpose()).transpose();
		SparseMatrix ret = new SparseMatrix(getRows(),matrix.getCols());
		for (int i = 0; i < sparseEntries; i++)
			for (int j = 0; j < ((SparseMatrix) matrix).sparseEntries; j++)
			{
				if(sparseXs[i] == ((SparseMatrix) matrix).sparseYs[i])
					ret.add(sparseValues[i]*((SparseMatrix) matrix).sparseValues[i],sparseYs[i],
						((SparseMatrix) matrix).sparseXs[i]);
			}
		return ret;
	}
	
	@Override
	public SparseMatrix getLowerTriangleMatrix()
	{
		SparseMatrix mat = new SparseMatrix(getRows(),getCols());
		for (int i = 0; i < sparseEntries; i++)
		{
			if(sparseYs[i] > sparseXs[i])
				mat.add(sparseValues[i],sparseYs[i],sparseXs[i]);
		}
		return mat;
	}
	
	@Override
	public SparseMatrix getUpperTriangleMatrix()
	{
		SparseMatrix mat = new SparseMatrix(getRows(),getCols());
		for (int i = 0; i < sparseEntries; i++)
		{
			if(sparseYs[i] < sparseXs[i])
				mat.add(sparseValues[i],sparseYs[i],sparseXs[i]);
		}
		return mat;
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
	@Override
	public SparseMatrix getDiagonalMatrix()
	{
		SparseMatrix mat = new SparseMatrix(getRows(),getCols());
		for (int i = 0; i < sparseEntries; i++)
		{
			if(sparseYs[i] == sparseXs[i])
				mat.add(sparseValues[i],sparseYs[i],sparseXs[i]);
		}
		return mat;
	}
}
