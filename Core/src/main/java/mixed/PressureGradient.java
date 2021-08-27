package mixed;

import linalg.CoordinateDenseMatrix;
import linalg.CoordinateVector;
import linalg.Tensor;

public class PressureGradient extends MixedGradient
{
	
	protected PressureGradient(final CoordinateVector pressureGradient)
	{
		super(pressureGradient.getLength());
		setPressureGradient(pressureGradient);
	}
	
	@Override
	public CoordinateDenseMatrix getVelocityGradient()
	{
		throw new IllegalStateException("Is not a velocity  gradient");
	}
	
	@Override
	public void addVelocityGradient(final CoordinateDenseMatrix velocityGradient)
	{
		throw new IllegalStateException("Is not a velocity  gradient");
	}
	
	@Override
	public double at(final int... coordinates)
	{
		if (coordinates[0] != 0) throw new IllegalStateException("Is not a velocity  gradient");
		return super.at(coordinates);
	}
	
	@Override
	public void set(final double value, final int... coordinates)
	{
		if (coordinates[0] != 0) throw new IllegalStateException("Is not a velocity  gradient");
		super.set(value, coordinates);
	}
	
	@Override
	public void add(final double value, final int... coordinates)
	{
		if (coordinates[0] != 0) throw new IllegalStateException("Is not a velocity  gradient");
		super.add(value, coordinates);
	}
	
	@Override
	public MixedGradient add(final Tensor other)
	{
		if (!(other instanceof MixedGradient))
			throw new IllegalArgumentException("can only add other mixedgradient");
		if (!getShape().equals(other.getShape())) throw new IllegalArgumentException("Wrong domain dimension");
		if (other instanceof PressureGradient) return new PressureGradient(
			getPressureGradient().add(((PressureGradient) other).getPressureGradient()));
		else if (other instanceof VelocityGradient)
		{
			final MixedGradient ret = new MixedGradient(getDomainDimension());
			ret.setPressureGradient(getPressureGradient());
			ret.setVelocityGradient(((VelocityGradient) other).getVelocityGradient());
			return ret;
		} else
		{
			final MixedGradient ret = new MixedGradient((MixedGradient) other);
			ret.addPressureGradient(getPressureGradient());
			return ret;
		}
	}
	
	@Override
	public PressureGradient mul(final double scalar)
	{
		return new PressureGradient(getPressureGradient().mul(scalar));
	}
}
