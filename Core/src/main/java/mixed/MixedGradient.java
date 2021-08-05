package mixed;

import basic.PerformanceArguments;
import linalg.*;

import java.util.List;

public class MixedGradient implements MutableMatrix
{
	private CoordinateVector pressureGradient;
	private CoordinateMatrix velocityGradient;
	
	protected MixedGradient(int domainDimension)
	{
		pressureGradient = new CoordinateVector(domainDimension);
		velocityGradient = new CoordinateMatrix(domainDimension, domainDimension);
	}
	
	public MixedGradient(MixedGradient mixedGradient)
	{
		this.pressureGradient = new CoordinateVector(mixedGradient.pressureGradient);
		this.velocityGradient = new CoordinateMatrix(mixedGradient.velocityGradient);
	}
	
	public CoordinateVector getPressureGradient()
	{
		return pressureGradient;
	}
	
	public void setPressureGradient(CoordinateVector pressureGradient)
	{
		this.pressureGradient = pressureGradient;
	}
	
	public void addPressureGradient(CoordinateVector pressureGradient)
	{
		this.setPressureGradient(getPressureGradient().add(pressureGradient));
	}
	
	public CoordinateMatrix getVelocityGradient()
	{
		return velocityGradient;
	}
	
	public void setVelocityGradient(CoordinateMatrix velocityGradient)
	{
		this.velocityGradient = velocityGradient;
	}
	
	public void addVelocityGradient(CoordinateMatrix velocityGradient)
	{
		this.setVelocityGradient(getVelocityGradient().add(velocityGradient));
	}
	
	@Override
	public double at(int... coordinates)
	{
		if (PerformanceArguments.getInstance().executeChecks) if (coordinates.length != 2)
			throw new IllegalArgumentException("Wrong number of coordinates");
		if (coordinates[0] == 0)
			return pressureGradient.at(coordinates[1]);
		else
			return velocityGradient.at(coordinates[0] - 1, coordinates[1]);
		
	}
	
	@Override
	public void set(double value, int... coordinates)
	{
		if (coordinates.length != 2)
			throw new IllegalArgumentException("Wrong number of coordinates");
		if (coordinates[0] == 0)
			pressureGradient.set(value, coordinates[1]);
		else
			velocityGradient.set(value, coordinates[0] - 1, coordinates[1]);
	}
	
	@Override
	public void add(double value, int... coordinates)
	{
		if (PerformanceArguments.getInstance().executeChecks) if (coordinates.length != 2)
			throw new IllegalArgumentException("Wrong number of coordinates");
		
		if (coordinates[0] == 0)
			pressureGradient.add(value, coordinates[1]);
		else
			velocityGradient.add(value, coordinates[0] - 1, coordinates[1]);
	}
	
	@Override
	public void addInPlace(Tensor other)
	{
	
	}
	
	@Override
	public void subInPlace(Tensor other)
	{
		MutableMatrix.super.subInPlace(other);
	}
	
	@Override
	public void mulInPlace(double scalar)
	{
	
	}
	
	
	@Override
	public MixedGradient add(Tensor other)
	{
		if (PerformanceArguments.getInstance().executeChecks)
		{
			if (!(other instanceof MixedGradient))
				throw new IllegalArgumentException("can only add other mixedgradient");
			if (!getShape().equals(other.getShape()))
				throw new IllegalArgumentException("Wrong domain dimension");
		}
		MixedGradient ret = new MixedGradient(this);
		if (!(other instanceof PressureGradient))
			ret.addVelocityGradient(((MixedGradient) other).getVelocityGradient());
		if (!(other instanceof VelocityGradient))
			ret.addPressureGradient(((MixedGradient) other).getPressureGradient());
		return ret;
	}
	
	public int getDomainDimension()
	{
		return pressureGradient.getLength();
	}
	
	@Override
	public MixedGradient mul(double scalar)
	{
		MixedGradient ret = new MixedGradient(getDomainDimension());
		ret.setVelocityGradient(getVelocityGradient().mul(scalar));
		ret.setPressureGradient(getPressureGradient().mul(scalar));
		return ret;
	}
	
	@Override
	public List<Vector> unfoldDimension(int dimension)
	{
		throw new UnsupportedOperationException("not implemented yet");
	}
	
	@Override
	public int getSparseEntryCount()
	{
		return 0;
	}
	
	@Override
	public boolean isSparse()
	{
		return false;
	}
	
	@Override
	public IntCoordinates getShape()
	{
		return new IntCoordinates(velocityGradient.getShape().get(0) + 1, velocityGradient.getShape().get(1) + 1);
	}
	
	@Override
	public Matrix transpose()
	{
		throw new UnsupportedOperationException("not implemented yet");
	}
	
	@Override
	public Vector mvMul(Vector vector)
	{
		throw new UnsupportedOperationException("not implemented yet");
	}
	
	@Override
	public Matrix mmMul(Matrix matrix)
	{
		throw new UnsupportedOperationException("not implemented yet");
	}
	
	@Override
	public void deleteRow(int row)
	{
		if(row == 0)
			pressureGradient.mulInPlace(0);
		else
			velocityGradient.deleteRow(row-1);
	}
	
	@Override
	public void deleteColumn(int column)
	{
		pressureGradient.set(0,column);
		velocityGradient.deleteColumn(column);
	}
}
