package linalg;

public class AffineTransformation
{
	public final CoordinateMatrix matrix;
	public final CoordinateVector vector;
	
	public AffineTransformation(CoordinateMatrix matrix, CoordinateVector vector)
	{
		this.matrix = matrix;
		this.vector = vector;
	}
	public CoordinateVector apply(CoordinateVector x)
	{
		return matrix.mvMul(x).add(vector);
	}
	
	@Override
	public String toString()
	{
		return "AffineTransformation{" +
			"matrix=" + matrix +
			", vector=" + vector +
			'}';
	}
}