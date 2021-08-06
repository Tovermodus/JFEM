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
	
	default VectorMultiplyable transpose()
	{
		VectorMultiplyable me = this;
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
			public Vector mvMul(Vector vector)
			{
				return me.tvMul(vector);
			}
			
			@Override
			public Vector tvMul(Vector vector)
			{
				return me.mvMul(vector);
			}
		};
	}
	
	default double powerIterationSymmetric()
	{
		if(PerformanceArguments.getInstance().executeChecks)
		{
			if(getVectorSize() != getTVectorSize())
				throw new IllegalStateException("Object is not square");
			if(IntStream.range(0, getVectorSize()).parallel()
				.mapToObj(i -> DenseVector.getUnitVector(getVectorSize(), i))
				.anyMatch(v -> !mvMul(v).almostEqual(tvMul(v))))
				throw new IllegalStateException("Not Symmetric");
		}
		DenseVector d;
		d = new DenseVector(getVectorSize());
		int component = 0;
		d.set(1,component++);
		while (this.mvMul(d).absMaxElement() < PerformanceArguments.getInstance().doubleTolerance*getVectorSize())
		{
			d = new DenseVector(getVectorSize());
			d.set(1,component++);
		}
		Vector b = d;
		Vector Ab;
		double old_lamb = 7.0;
		double new_lamb = 6.0;
		while (!DoubleCompare.almostEqualAfterOps(old_lamb/new_lamb, 1, this.getVectorSize()))
		{
			old_lamb = new_lamb;
			Ab = this.mvMul(b);
			new_lamb = b.inner(Ab)/Math.pow(b.euclidianNorm(),2);
			b = Ab.mul(1./Ab.euclidianNorm());
		}
		return new_lamb;
	}
}
