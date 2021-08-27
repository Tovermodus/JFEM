package mixed;

import basic.PerformanceArguments;
import linalg.*;

import java.util.List;

public class MixedGradient implements MutableMatrix
{
	private CoordinateVector pressureGradient;
	private CoordinateDenseMatrix velocityGradient;
	
	protected MixedGradient(final int domainDimension)
	{
		pressureGradient = new CoordinateVector(domainDimension);
		velocityGradient = new CoordinateDenseMatrix(domainDimension, domainDimension);
	}
	
	public MixedGradient(final MixedGradient mixedGradient)
	{
		this.pressureGradient = new CoordinateVector(mixedGradient.pressureGradient);
		this.velocityGradient = new CoordinateDenseMatrix(mixedGradient.velocityGradient);
	}
	
	public CoordinateVector getPressureGradient()
	{
		return pressureGradient;
	}
	
	public void setPressureGradient(final CoordinateVector pressureGradient)
	{
		this.pressureGradient = pressureGradient;
	}
	
	public void addPressureGradient(final CoordinateVector pressureGradient)
	{
		this.setPressureGradient(getPressureGradient().add(pressureGradient));
	}
	
	public CoordinateDenseMatrix getVelocityGradient()
	{
		return velocityGradient;
	}
	
	public void setVelocityGradient(final CoordinateMatrix velocityGradient)
	{
		if (velocityGradient instanceof CoordinateDenseMatrix)
			this.velocityGradient = (CoordinateDenseMatrix) velocityGradient;
		else this.velocityGradient = new CoordinateDenseMatrix(velocityGradient);
	}
	
	public void addVelocityGradient(final CoordinateDenseMatrix velocityGradient)
	{
		this.setVelocityGradient(getVelocityGradient().add(velocityGradient));
	}
	
	@Override
	public double at(final int... coordinates)
	{
		if (PerformanceArguments.getInstance().executeChecks)
			if (coordinates.length != 2) throw new IllegalArgumentException("Wrong number of coordinates");
		if (coordinates[0] == 0) return pressureGradient.at(coordinates[1]);
		else return velocityGradient.at(coordinates[0] - 1, coordinates[1]);
	}
	
	@Override
	public void set(final double value, final int... coordinates)
	{
		if (coordinates.length != 2) throw new IllegalArgumentException("Wrong number of coordinates");
		if (coordinates[0] == 0) pressureGradient.set(value, coordinates[1]);
		else velocityGradient.set(value, coordinates[0] - 1, coordinates[1]);
	}
	
	@Override
	public void add(final double value, final int... coordinates)
	{
		if (PerformanceArguments.getInstance().executeChecks)
			if (coordinates.length != 2) throw new IllegalArgumentException("Wrong number of coordinates");
		
		if (coordinates[0] == 0) pressureGradient.add(value, coordinates[1]);
		else velocityGradient.add(value, coordinates[0] - 1, coordinates[1]);
	}
	
	@Override
	public void addInPlace(final Tensor other)
	{
	
	}
	
	@Override
	public void subInPlace(final Tensor other)
	{
		MutableMatrix.super.subInPlace(other);
	}
	
	@Override
	public void mulInPlace(final double scalar)
	{
	
	}
	
	@Override
	public MixedGradient add(final Tensor other)
	{
		if (PerformanceArguments.getInstance().executeChecks)
		{
			if (!(other instanceof MixedGradient))
				throw new IllegalArgumentException("can only add other mixedgradient");
			if (!getShape().equals(other.getShape()))
				throw new IllegalArgumentException("Wrong domain dimension");
		}
		final MixedGradient ret = new MixedGradient(this);
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
	public MixedGradient mul(final double scalar)
	{
		final MixedGradient ret = new MixedGradient(getDomainDimension());
		ret.setVelocityGradient(getVelocityGradient().mul(scalar));
		ret.setPressureGradient(getPressureGradient().mul(scalar));
		return ret;
	}
	
	@Override
	public List<Vector> unfoldDimension(final int dimension)
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
		return new IntCoordinates(velocityGradient.getShape().get(0) + 1,
		                          velocityGradient.getShape().get(1) + 1);
	}
	
	@Override
	public Matrix transpose()
	{
		throw new UnsupportedOperationException("not implemented yet");
	}
	
	@Override
	public Vector mvMul(final Vector vector)
	{
		throw new UnsupportedOperationException("not implemented yet");
	}
	
	@Override
	public Matrix mmMul(final Matrix matrix)
	{
		throw new UnsupportedOperationException("not implemented yet");
	}
	
	@Override
	public void deleteRow(final int row)
	{
		if (row == 0) pressureGradient.mulInPlace(0);
		else velocityGradient.deleteRow(row - 1);
	}
	
	@Override
	public void deleteColumn(final int column)
	{
		pressureGradient.set(0, column);
		velocityGradient.deleteColumn(column);
	}
}
