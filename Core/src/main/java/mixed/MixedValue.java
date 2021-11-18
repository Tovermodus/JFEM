package mixed;

import basic.PerformanceArguments;
import linalg.*;

public class MixedValue
	implements MutableVector
{
	private double pressure;
	private CoordinateVector velocity;
	
	public MixedValue(final MixedValue mv)
	{
		pressure = mv.getPressure();
		velocity = new CoordinateVector(mv.getVelocity());
	}
	
	public MixedValue(final double pressure, final CoordinateVector velocity)
	{
		this.pressure = pressure;
		this.velocity = velocity;
	}
	
	protected MixedValue(final int domainDimension)
	{
		pressure = 0;
		velocity = new CoordinateVector(domainDimension);
	}
	
	@Override
	public double at(final int... coordinates)
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
	public void set(final double value, final int... coordinates)
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
	public void add(final double value, final int... coordinates)
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
	public void addInPlace(final Tensor other)
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
	public void mulInPlace(final double scalar)
	{
		pressure *= scalar;
		velocity.mulInPlace(scalar);
	}
	
	@Override
	public MixedValue add(final Tensor other)
	{
		if (PerformanceArguments.getInstance().executeChecks)
		{
			if (!(other instanceof MixedValue))
				throw new IllegalArgumentException("Cant add MixedValueVector to different Vector");
			if (!getShape().equals(other.getShape()))
				throw new IllegalArgumentException("Vectors are of different size");
		}
		final MixedValue ret = new MixedValue(this);
		if (!(other instanceof PressureValue))
			ret.addVelocity(((MixedValue) other).getVelocity());
		if (!(other instanceof VelocityValue))
			ret.addPressure(((MixedValue) other).getPressure());
		return ret;
	}
	
	@Override
	public MixedValue sub(final Tensor other)
	{
		return this.add(other.mul(-1));
	}
	
	@Override
	public MixedValue mul(final double scalar)
	{
		final MixedValue ret = new MixedValue(getDomainDimension());
		ret.setPressure(getPressure() * scalar);
		ret.setVelocity(getVelocity().mul(scalar));
		return ret;
	}
	
	@Override
	public Matrix outer(final Vector other)
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
	
	public void setPressure(final double pressure)
	{
		this.pressure = pressure;
	}
	
	public void setVelocity(final CoordinateVector velocity)
	{
		this.velocity = new CoordinateVector(velocity);
	}
	
	public void addPressure(final double pressure)
	{
		this.pressure += pressure;
	}
	
	public void addVelocity(final CoordinateVector velocity)
	{
		setVelocity(getVelocity().add(velocity));
	}
	
	@Override
	public String toString()
	{
		return "Mixed, P: " + getPressure() + ", V: " + getVelocity().toString();
	}
}
