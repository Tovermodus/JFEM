package linalg;

public class CoordinateMatrix extends DenseMatrix
{
	public CoordinateMatrix(int i, int j)
	{
		super(i, j);
		if(i > 3|| j>3)
			throw new IllegalArgumentException("only 1D, 2D and 3D supported");
	
	}
	
	public CoordinateMatrix(Matrix matrix)
	{
		super(matrix);
		if(matrix.getShape().get(0) > 3|| matrix.getShape().get(1) > 3)
			throw new IllegalArgumentException("only 1D, 2D and 3D supported");
	}
	
	@Override
	public CoordinateVector mvMul(Vector vector)
	{
		return (CoordinateVector) super.mvMul(vector);
	}
	
	@Override
	public CoordinateMatrix add(Tensor other)
	{
		if(!getShape().equals(other.getShape()))
			throw new IllegalArgumentException("Incompatible sizes");
		CoordinateMatrix ret = new CoordinateMatrix(this);
		for(int i = 0; i < ret.getRows(); i++)
			for(int j = 0; j < ret.getCols(); j++)
				ret.add(other.at(i,j),i,j);
		return ret;
	}
	
	@Override
	public CoordinateMatrix sub(Tensor other){
		
		if(!getShape().equals(other.getShape()))
			throw new IllegalArgumentException("Incompatible sizes");
		CoordinateMatrix ret = new CoordinateMatrix(this);
		for(int i = 0; i < ret.getRows(); i++)
			for(int j = 0; j < ret.getCols(); j++)
				ret.add(-other.at(i,j),i,j);
		return ret;
	}
	
	@Override
	public CoordinateMatrix mul(double scalar)
	{
		
		CoordinateMatrix ret = new CoordinateMatrix(entries.length, entries[0].length);
		for(int i = 0; i < ret.getRows(); i++)
			for(int j = 0; j < ret.getCols(); j++)
				ret.set(scalar*at(i,j),i,j);
		return ret;
	}
}
