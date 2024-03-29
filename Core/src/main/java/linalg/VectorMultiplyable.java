package linalg;

import basic.DoubleCompare;
import basic.PerformanceArguments;

import java.util.stream.IntStream;

public interface VectorMultiplyable
{
	int getVectorSize();
	
	int getTVectorSize();
	
	Vector mvMul(Vector vector);
	
	Vector tvMul(Vector vector);
	
	static VectorMultiplyable concatenate(final VectorMultiplyable v1, final VectorMultiplyable v2)
	{
		if (v1.getVectorSize() != v2.getTVectorSize())
			throw new IllegalArgumentException("Sizes dont match" + v1.getVectorSize() + "!=" + v2.getTVectorSize());
		return new VectorMultiplyable()
		{
			@Override
			public int getVectorSize()
			{
				return v2.getVectorSize();
			}
			
			@Override
			public int getTVectorSize()
			{
				return v1.getTVectorSize();
			}
			
			@Override
			public Vector mvMul(final Vector vector)
			{
				return v1.mvMul(v2.mvMul(vector));
			}
			
			@Override
			public Vector tvMul(final Vector vector)
			{
				return v2.tvMul(v1.tvMul(vector));
			}
		};
	}
	
	static VectorMultiplyable concatenateTranspose(final VectorMultiplyable v1, final VectorMultiplyable v2)
	{
		if (v1.getVectorSize() != v2.getVectorSize())
			throw new IllegalArgumentException("Sizes dont match" + v1.getVectorSize() + "!=" + v2.getVectorSize());
		return new VectorMultiplyable()
		{
			@Override
			public int getVectorSize()
			{
				return v2.getTVectorSize();
			}
			
			@Override
			public int getTVectorSize()
			{
				return v1.getTVectorSize();
			}
			
			@Override
			public Vector mvMul(final Vector vector)
			{
				return v1.mvMul(v2.tvMul(vector));
			}
			
			@Override
			public Vector tvMul(final Vector vector)
			{
				return v2.mvMul(v1.tvMul(vector));
			}
		};
	}
	
	static VectorMultiplyable transposeConcatenate(final VectorMultiplyable v1, final VectorMultiplyable v2)
	{
		if (v1.getTVectorSize() != v2.getTVectorSize())
			throw new IllegalArgumentException("Sizes dont match");
		return new VectorMultiplyable()
		{
			@Override
			public int getVectorSize()
			{
				return v2.getVectorSize();
			}
			
			@Override
			public int getTVectorSize()
			{
				return v1.getVectorSize();
			}
			
			@Override
			public Vector mvMul(final Vector vector)
			{
				return v1.tvMul(v2.mvMul(vector));
			}
			
			@Override
			public Vector tvMul(final Vector vector)
			{
				return v2.tvMul(v1.mvMul(vector));
			}
		};
	}
	
	default VectorMultiplyable addVm(final VectorMultiplyable other)
	{
		if (other.getVectorSize() != getVectorSize())
			throw new IllegalArgumentException("Sizes dont match");
		if (getTVectorSize() != other.getTVectorSize())
			throw new IllegalArgumentException("Sizes dont match");
		assert other.getTVectorSize() == getTVectorSize();
		assert other.getVectorSize() == getVectorSize();
		final VectorMultiplyable me = this;
		return new VectorMultiplyable()
		{
			@Override
			public int getVectorSize()
			{
				return me.getTVectorSize();
			}
			
			@Override
			public int getTVectorSize()
			{
				return me.getVectorSize();
			}
			
			@Override
			public Vector mvMul(final Vector vector)
			{
				return me.mvMul(vector)
				         .add(other.mvMul(vector));
			}
			
			@Override
			public Vector tvMul(final Vector vector)
			{
				return me.tvMul(vector)
				         .add(other.tvMul(vector));
			}
		};
	}
	
	default VectorMultiplyable mulVm(final double scalar)
	{
		final VectorMultiplyable me = this;
		return new VectorMultiplyable()
		{
			@Override
			public int getVectorSize()
			{
				return me.getTVectorSize();
			}
			
			@Override
			public int getTVectorSize()
			{
				return me.getVectorSize();
			}
			
			@Override
			public Vector mvMul(final Vector vector)
			{
				return me.mvMul(vector)
				         .mul(scalar);
			}
			
			@Override
			public Vector tvMul(final Vector vector)
			{
				return me.tvMul(vector)
				         .mul(scalar);
			}
		};
	}
	
	default VectorMultiplyable transpose()
	{
		final VectorMultiplyable me = this;
		return new VectorMultiplyable()
		{
			@Override
			public int getVectorSize()
			{
				return me.getTVectorSize();
			}
			
			@Override
			public int getTVectorSize()
			{
				return me.getVectorSize();
			}
			
			@Override
			public Vector mvMul(final Vector vector)
			{
				return me.tvMul(vector);
			}
			
			@Override
			public Vector tvMul(final Vector vector)
			{
				return me.mvMul(vector);
			}
		};
	}
	
	default double powerIterationSymmetric()
	{
		if (PerformanceArguments.getInstance().executeChecks)
		{
			if (getVectorSize() != getTVectorSize())
				throw new IllegalStateException("Object is not square");
			if (IntStream.range(0, getVectorSize())
			             .parallel()
			             .mapToObj(i -> DenseVector.getUnitVector(getVectorSize(), i))
			             .anyMatch(v -> !mvMul(v).almostEqual(tvMul(v))))
				throw new IllegalStateException("Not Symmetric");
		}
		DenseVector d;
		d = new DenseVector(getVectorSize());
		for (int i = 0; i < d.size(); i++)
			d.set(Math.random(), i);
		int component = 0;
		d.add(1, component++);
		while (this.mvMul(d)
		           .absMaxElement() < PerformanceArguments.getInstance().doubleTolerance * getVectorSize())
		{
			d = new DenseVector(getVectorSize());
			d.set(1, component++);
		}
		d = d.mul(1. / d.euclidianNorm());
		Vector b = d;
		Vector Ab;
		double old_lamb = 7.0;
		double new_lamb = 6.0;
		while (!DoubleCompare.almostEqualAfterOps(old_lamb / new_lamb, 1, this.getVectorSize()))
		{
			old_lamb = new_lamb;
			Ab = this.mvMul(b);
			new_lamb = b.inner(Ab) / Math.pow(b.euclidianNorm(), 2);
			b = Ab.mul(1. / Ab.euclidianNorm());
		}
		return new_lamb;
	}
	
	default double powerIterationNonSymmetric()
	{
		if (PerformanceArguments.getInstance().executeChecks)
		{
			if (getVectorSize() != getTVectorSize())
				throw new IllegalStateException("Object is not square");
		}
		DenseVector d;
		d = new DenseVector(getVectorSize());
		for (int i = 0; i < d.size(); i++)
			d.set(Math.random(), i);
		int component = 0;
		d.add(1, component++);
		while (this.mvMul(d)
		           .absMaxElement() < PerformanceArguments.getInstance().doubleTolerance * getVectorSize())
		{
			d = new DenseVector(getVectorSize());
			d.set(1, component++);
		}
		d = d.mul(1. / d.euclidianNorm());
		Vector b = d;
		Vector Ab;
		double old_lamb = 7.0;
		double new_lamb = 6.0;
		while (!DoubleCompare.almostEqualAfterOps(old_lamb / new_lamb, 1, this.getVectorSize()))
		{
			old_lamb = new_lamb;
			Ab = this.tvMul(this.mvMul(b));
			new_lamb = b.inner(Ab) / Math.pow(b.euclidianNorm(), 2);
			b = Ab.mul(1. / Ab.euclidianNorm());
		}
		return Math.sqrt(new_lamb);
	}
}
