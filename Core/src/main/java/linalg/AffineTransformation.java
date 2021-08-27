package linalg;

public class AffineTransformation
{
	public final CoordinateDenseMatrix matrix;
	public final CoordinateVector vector;
	
	public AffineTransformation(final CoordinateDenseMatrix matrix, final CoordinateVector vector)
	{
		this.matrix = matrix;
		this.vector = vector;
	}
	
	public CoordinateVector apply(final CoordinateVector x)
	{
		return matrix.mvMul(x).add(vector);
	}
	
	public CoordinateVector applyInverse(final CoordinateVector x)
	{
		return matrix.solve(x.sub(vector));
	}
	
	@Override
	public String toString()
	{
		return "AffineTransformation{" + "matrix=" + matrix + ", vector=" + vector + '}';
	}
}
