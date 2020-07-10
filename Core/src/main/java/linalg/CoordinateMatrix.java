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
	
}
