package linalg;

import basic.PerformanceArguments;

import java.util.List;
import java.util.Map;
import java.util.TreeMap;

public interface Matrix extends Tensor, VectorMultiplyable
{
	@Override
	Matrix add(Tensor other);
	
	@Override
	default Matrix sub(Tensor other)
	{
		if (!getShape().equals(other.getShape()))
			throw new IllegalArgumentException("Matrices are of different size");
		return add(other.mul(-1.));
	}
	
	
	@Override
	Matrix mul(double scalar);
	
	default int getRows()
	{
		return getShape().get(0);
	}
	default int getCols()
	{
		return getShape().get(1);
	}
	@Override
	List<Vector> unfoldDimension(int dimension);
	
	Matrix transpose();
	
	@Override
	Vector mvMul(Vector vector);
	
	default Vector tvMul(Vector vector)
	{
		if (PerformanceArguments.getInstance().executeChecks)
			if (getRows() != vector.getLength())
				throw new IllegalArgumentException("Incompatible sizes");
		return transpose().mvMul(vector);
	}
	
	default double frobeniusInner(Matrix other)
	{
		throw new UnsupportedOperationException("not yet implemented");
	}
	
	Matrix mmMul(Matrix matrix);
	
	default Matrix tmMul(Matrix matrix)
	{
		if (PerformanceArguments.getInstance().executeChecks)
			if (getRows() != (matrix.getRows()))
				throw new IllegalArgumentException("Incompatible sizes");
		return transpose().mmMul(matrix);
	}
	
	default Matrix mtMul(Matrix matrix)
	{
		if (PerformanceArguments.getInstance().executeChecks)
			if (getCols() != (matrix.getCols()))
				throw new IllegalArgumentException("Incompatible sizes");
		return mmMul(matrix.transpose());
	}
	
	@Override
	default String printFormatted(double... tol)
	{
		if (PerformanceArguments.getInstance().executeChecks)
			if (tol.length > 1)
				throw new IllegalArgumentException("Only one Tolerance accepted");
		if (tol.length == 0)
			tol = new double[]{1e-14};
		String ret = "[";
		for(int i = 0; i < getRows(); i++)
		{
			for(int j = 0; j < getCols(); j++)
			{
				if (Math.abs(at(i,j)) < tol[0])
					ret = ret.concat("<>........  ");
				else
				{
					if (at(i,j) >= 0)
						ret = ret.concat("+");
					ret = ret.concat(String.format("%6.3e", at(i,j)) + "  ");
				}
			}
			ret = ret.concat("\n ");
		}
		ret = ret.concat("]");
		return ret;
	}
}
