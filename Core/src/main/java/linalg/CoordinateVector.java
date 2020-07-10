package linalg;

public class CoordinateVector extends DenseVector
{
	public CoordinateVector(Vector v)
	{
		super(v);
		if(v.getLength()>3)
			throw new IllegalArgumentException("only 1D, 2D and 3D supported");
	}
	public CoordinateVector(int d)
	{
		super(d);
		if(d > 3)
			throw new IllegalArgumentException("only 1D, 2D and 3D supported");
	}
	public static CoordinateVector getE(int d, int index)
	{
		CoordinateVector ret = new CoordinateVector(d);
		ret.set(index, 1);
		return ret;
	}
	public static CoordinateVector fromValues(double ... values)
	{
		CoordinateVector ret = new CoordinateVector(values.length);
		for (int i = 0; i < values.length; i++)
		{
			ret.set(values[i],i);
		}
		return ret;
	}
	public double x()
	{
		return entries[0];
	}
	public double y()
	{
		return entries[1];
	}
	public double z()
	{
		return entries[2];
	}
}
