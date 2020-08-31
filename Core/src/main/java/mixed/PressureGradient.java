package mixed;

import linalg.CoordinateMatrix;
import linalg.CoordinateVector;
import linalg.Tensor;

public class PressureGradient extends MixedGradient
{
	
	protected PressureGradient(CoordinateVector pressureGradient)
	{
		super(pressureGradient.getLength());
		setPressureGradient(pressureGradient);
	}
	
	@Override
	public CoordinateMatrix getVelocityGradient()
	{
		throw new IllegalStateException("Is not a velocity  gradient");
	}
	
	@Override
	public void setVelocityGradient(CoordinateMatrix velocityGradient)
	{
		throw new IllegalStateException("Is not a velocity  gradient");
	}
	
	@Override
	public void addVelocityGradient(CoordinateMatrix velocityGradient)
	{
		throw new IllegalStateException("Is not a velocity  gradient");
	}
	
	@Override
	public double at(int... coordinates)
	{
		if(coordinates[0] != 0)
			throw new IllegalStateException("Is not a velocity  gradient");
		return super.at(coordinates);
	}
	
	@Override
	public void set(double value, int... coordinates)
	{
		if(coordinates[0] != 0)
			throw new IllegalStateException("Is not a velocity  gradient");
		super.set(value, coordinates);
	}
	
	@Override
	public void add(double value, int... coordinates)
	{
		if(coordinates[0] != 0)
			throw new IllegalStateException("Is not a velocity  gradient");
		super.add(value, coordinates);
	}
	
	@Override
	public MixedGradient add(Tensor other)
	{
		if(!(other instanceof MixedGradient))
			throw new IllegalArgumentException("can only add other mixedgradient");
		if (!getShape().equals(other.getShape()))
			throw new IllegalArgumentException("Wrong domain dimension");
		if(other instanceof PressureGradient)
			return new PressureGradient(getPressureGradient().add(((PressureGradient) other).getPressureGradient()));
		else if(other instanceof VelocityGradient)
		{
			MixedGradient ret = new MixedGradient(getDomainDimension());
			ret.setPressureGradient(getPressureGradient());
			ret.setVelocityGradient(((VelocityGradient) other).getVelocityGradient());
			return ret;
		}
		else
		{
			MixedGradient ret = new MixedGradient((MixedGradient) other);
			ret.addPressureGradient(getPressureGradient());
			return ret;
		}
	}
	
	@Override
	public PressureGradient mul(double scalar)
	{
		return new PressureGradient(getPressureGradient().mul(scalar));
	}
}
