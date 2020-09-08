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
	public static CoordinateVector getUnitVector(int d, int index)
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
	
	
	public CoordinateMatrix outer(CoordinateVector other)
	{
		return new CoordinateMatrix(super.outer(other));
	}
	
	@Override
	public CoordinateVector add(Tensor other)
	{
		if(!(other instanceof CoordinateVector))
			throw new IllegalArgumentException("can't add coordinates and other vectors");
		if(!getShape().equals(other.getShape()))
			throw new IllegalArgumentException("Vectors are of different size");
		CoordinateVector ret = new CoordinateVector(this);
		for (int i = 0; i < getLength(); i++)
			ret.entries[i] = entries[i]+((DenseVector) other).entries[i];
		return ret;
	}
	@Override
	public CoordinateVector sub(Tensor other)
	{
		if(!(other instanceof CoordinateVector))
			throw new IllegalArgumentException("can't add coordinates and other vectors");
		if(!getShape().equals(other.getShape()))
			throw new IllegalArgumentException("Vectors are of different size");
		CoordinateVector ret = new CoordinateVector(this);
		for (int i = 0; i < getLength(); i++)
			ret.entries[i] = entries[i]-((DenseVector) other).entries[i];
		return ret;
	}
	
	@Override
	public CoordinateVector mul(double scalar)
	{
		CoordinateVector ret = new CoordinateVector(entries.length);
		for (int i = 0; i < getLength(); i++)
		{
			ret.set(at(i)*scalar,i);
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
