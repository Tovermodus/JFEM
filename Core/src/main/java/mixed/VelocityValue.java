package mixed;

import linalg.CoordinateVector;
import linalg.Tensor;

public class VelocityValue
	extends MixedValue
{
	public VelocityValue(final CoordinateVector velocity)
	{
		super(0, new CoordinateVector(velocity));
	}
	
	private VelocityValue(final CoordinateVector velocity, final boolean wrap)
	{
		super(0, velocity);
	}
	
	public VelocityValue(final int domainDimension)
	{
		super(domainDimension);
	}
	
	public static VelocityValue wrap(final CoordinateVector c)
	{
		return new VelocityValue(c, true);
	}
	
	@Override
	public double at(final int... coordinates)
	{
		if (coordinates.length != 1)
			throw new IllegalArgumentException("Wrong number of coordinates");
		if (coordinates[0] == 0)
			throw new IllegalStateException("Not a pressure vector");
		return getVelocity().at(coordinates[0] - 1);
	}
	
	@Override
	public void set(final double value, final int... coordinates)
	{
		if (coordinates.length != 1)
			throw new IllegalArgumentException("Wrong number of coordinates");
		if (coordinates[0] == 0)
			throw new IllegalStateException("Not a pressure vector");
		getVelocity().set(value, coordinates[0] - 1);
	}
	
	@Override
	public void add(final double value, final int... coordinates)
	{
		if (coordinates.length != 1)
			throw new IllegalArgumentException("Wrong number of coordinates");
		if (coordinates[0] == 0)
			throw new IllegalStateException("Not a pressure vector");
		getVelocity().add(value, coordinates[0] - 1);
	}
	
	@Override
	public MixedValue add(final Tensor other)
	{
		if (other instanceof VelocityValue)
		{
			if (getDomainDimension() != ((VelocityValue) other).getDomainDimension())
				throw new IllegalArgumentException("Wrong domain dimension");
			return new VelocityValue(getVelocity().add(((VelocityValue) other).getVelocity()));
		} else if (other instanceof PressureValue)
		{
			final MixedValue ret = new MixedValue(getDomainDimension());
			ret.setPressure(((PressureValue) other).getPressure());
			ret.setVelocity(getVelocity());
			return ret;
		} else if (other instanceof MixedValue)
		{
			if (getDomainDimension() != ((MixedValue) other).getDomainDimension())
				throw new IllegalArgumentException("Wrong domain dimension");
			final MixedValue ret = new MixedValue((MixedValue) other);
			ret.addVelocity(getVelocity());
			return ret;
		} else
			return super.add(other);
	}
	
	@Override
	public VelocityValue mul(final double scalar)
	{
		return new VelocityValue(getVelocity().mul(scalar));
	}
	
	@Override
	public double getPressure()
	{
		throw new IllegalStateException("not a pressure vector");
	}
	
	@Override
	public void setPressure(final double pressure)
	{
		throw new IllegalStateException("not a pressure vector");
	}
	
	@Override
	public void addPressure(final double pressure)
	{
		throw new IllegalStateException("not a pressure vector");
	}
	
	@Override
	public String toString()
	{
		return "Velocity: " + getVelocity().toString();
	}
}
