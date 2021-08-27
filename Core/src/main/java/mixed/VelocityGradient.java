package mixed;

import linalg.CoordinateMatrix;
import linalg.CoordinateVector;
import linalg.Tensor;

public class VelocityGradient extends MixedGradient
{
	public VelocityGradient(final CoordinateMatrix velocityGradient)
	{
		super(velocityGradient.getCols());
		setVelocityGradient(velocityGradient);
	}
	
	@Override
	public CoordinateVector getPressureGradient()
	{
		throw new IllegalStateException("Is not a pressure  gradient");
	}
	
	@Override
	public void setPressureGradient(final CoordinateVector pressureGradient)
	{
		throw new IllegalStateException("Is not a pressure  gradient");
	}
	
	@Override
	public void addPressureGradient(final CoordinateVector pressureGradient)
	{
		throw new IllegalStateException("Is not a pressure  gradient");
	}
	
	@Override
	public double at(final int... coordinates)
	{
		if (coordinates[0] == 0) throw new IllegalStateException("Is not a pressure  gradient");
		return super.at(coordinates);
	}
	
	@Override
	public void set(final double value, final int... coordinates)
	{
		if (coordinates[0] == 0) throw new IllegalStateException("Is not a pressure  gradient");
		super.set(value, coordinates);
	}
	
	@Override
	public void add(final double value, final int... coordinates)
	{
		if (coordinates[0] == 0) throw new IllegalStateException("Is not a pressure  gradient");
		super.add(value, coordinates);
	}
	
	@Override
	public MixedGradient add(final Tensor other)
	{
		if (!(other instanceof MixedGradient))
			throw new IllegalArgumentException("can only add other mixedgradient");
		if (!getShape().equals(other.getShape())) throw new IllegalArgumentException("Wrong domain dimension");
		if (other instanceof VelocityGradient) return new VelocityGradient(
			getVelocityGradient().add(((VelocityGradient) other).getVelocityGradient()));
		else if (other instanceof PressureGradient)
		{
			final MixedGradient ret = new MixedGradient(getDomainDimension());
			ret.setPressureGradient(((PressureGradient) other).getPressureGradient());
			ret.setVelocityGradient(getVelocityGradient());
			return ret;
		} else
		{
			final MixedGradient ret = new MixedGradient((MixedGradient) other);
			ret.addVelocityGradient(getVelocityGradient());
			return ret;
		}
	}
	
	@Override
	public VelocityGradient mul(final double scalar)
	{
		return new VelocityGradient(getVelocityGradient().mul(scalar));
	}
}
