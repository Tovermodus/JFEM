package mixed;

import basic.PerformanceArguments;
import linalg.*;


public class MixedValue implements MutableVector
{
	private double pressure;
	private CoordinateVector velocity;
	
	public MixedValue(MixedValue mv)
	{
		pressure = mv.getPressure();
		velocity = new CoordinateVector(mv.getVelocity());
	}
	
	protected MixedValue(int d)
	{
		pressure = 0;
		velocity = new CoordinateVector(d);
	}
	
	@Override
	public double at(int... coordinates)
	{
		
		if (PerformanceArguments.getInstance().executeChecks)
			if (coordinates.length != 1)
				throw new IllegalArgumentException("Wrong number of coordinates");
		if (coordinates[0] == 0)
			return getPressure();
		else
			return getVelocity().at(coordinates[0] - 1);
	}
	
	@Override
	public void set(double value, int... coordinates)
	{
		if (PerformanceArguments.getInstance().executeChecks)
			if (coordinates.length != 1)
				throw new IllegalArgumentException("Wrong number of coordinates");
		if (coordinates[0] == 0)
		{
			setPressure(value);
		} else
		{
			getVelocity().set(value, coordinates[0] - 1);
		}
	}
	
	@Override
	public void add(double value, int... coordinates)
	{
		if (PerformanceArguments.getInstance().executeChecks)
			if (coordinates.length != 1)
				throw new IllegalArgumentException("Wrong number of coordinates");
		if (coordinates[0] == 0)
		{
			setPressure(getPressure() + value);
			
		} else
		{
			getVelocity().add(value, coordinates[0] - 1);
		}
	}
	
	@Override
	public void addInPlace(Tensor other)
	{
		if (PerformanceArguments.getInstance().executeChecks)
		{
			if (!(other instanceof MixedValue))
				throw new IllegalArgumentException("Cant add MixedValueVector to different Vector");
			if (!getShape().equals(other.getShape()))
				throw new IllegalArgumentException("Vectors are of different size");
		}
		pressure += ((MixedValue) other).pressure;
		velocity.addInPlace(((MixedValue) other).velocity);
	}
	
	@Override
	public void mulInPlace(double scalar)
	{
		pressure *= scalar;
		velocity.mulInPlace(scalar);
	}
	
	@Override
	public MixedValue add(Tensor other)
	{
		if (PerformanceArguments.getInstance().executeChecks)
		{
			if (!(other instanceof MixedValue))
				throw new IllegalArgumentException("Cant add MixedValueVector to different Vector");
			if (!getShape().equals(other.getShape()))
				throw new IllegalArgumentException("Vectors are of different size");
		}
		MixedValue ret = new MixedValue(this);
		if (!(other instanceof PressureValue))
			ret.addVelocity(((MixedValue) other).getVelocity());
		if (!(other instanceof VelocityValue))
			ret.addPressure(((MixedValue) other).getPressure());
		return ret;
	}
	
	@Override
	public MixedValue mul(double scalar)
	{
		MixedValue ret = new MixedValue(getDomainDimension());
		ret.setPressure(getPressure() * scalar);
		ret.setVelocity(getVelocity().mul(scalar));
		return ret;
	}
	
	@Override
	public Matrix outer(Vector other)
	{
		throw new UnsupportedOperationException("outer makes no sense");
	}
	
	public CoordinateVector getVelocity()
	{
		return velocity;
	}
	
	public double getPressure()
	{
		return pressure;
	}
	
	
	@Override
	public int getSparseEntryCount()
	{
		return getLength();
	}
	
	@Override
	public boolean isSparse()
	{
		return false;
	}
	
	@Override
	public IntCoordinates getShape()
	{
		return new IntCoordinates(1 + velocity.getLength());
	}
	
	public IntCoordinates getVelocityShape()
	{
		return getVelocity().getShape();
	}
	
	public int getDomainDimension()
	{
		return getVelocity().getLength();
	}
	
	public void setPressure(double pressure)
	{
		this.pressure = pressure;
	}
	
	public void setVelocity(CoordinateVector velocity)
	{
		this.velocity = new CoordinateVector(velocity);
	}
	
	public void addPressure(double pressure)
	{
		this.pressure += pressure;
	}
	
	public void addVelocity(CoordinateVector velocity)
	{
		setVelocity(getVelocity().add(velocity));
	}
	
	@Override
	public String toString()
	{
		return "Mixed, P: " + getPressure() + ", V: " + getVelocity().toString();
	}
}
