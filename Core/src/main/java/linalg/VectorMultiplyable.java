package linalg;

import java.util.List;

public interface VectorMultiplyable
{
	int getVectorSize();
	Vector mvMul(Vector vector);
	
	default double powerIteration()
	{
		DenseVector d;
		d = new DenseVector(getVectorSize());
		int component = 0;
		d.set(1,component++);
		while (this.mvMul(d).absMaxElement() < 1e-14)
		{
			d = new DenseVector(getVectorSize());
			d.set(1,component++);
		}
		Vector b = d;
		Vector Ab;
		double old_lamb = 7.0;
		double new_lamb = 6.0;
		while (Math.abs(old_lamb / new_lamb - 1) > 1e-6)
		{
			old_lamb = new_lamb;
			Ab = this.mvMul(b);
			new_lamb = b.inner(Ab)/Math.pow(b.euclidianNorm(),2);
			b = Ab.mul(1./Ab.euclidianNorm());
		}
		return new_lamb;
	}	
}
