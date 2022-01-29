package linalg;

import java.util.ArrayList;
import java.util.List;

public class Arnoldi
{
	List<Vector> q;
	List<Vector> Av;
	List<Vector> v;
	final VectorMultiplyable A;
	
	public Arnoldi(final Vector initial, final VectorMultiplyable A)
	{
		this.A = A;
		q = new ArrayList<>();
		
		Av = new ArrayList<>();
		v = new ArrayList<>();
		
		q.add(initial);
		final Vector w0 = latestQNormalized();
		addQ();
		final Vector Aw0 = latestQ();
		final Vector v0 = w0.mul(1. / Aw0.euclidianNorm());
		final Vector Av0 = Aw0.mul(1. / Aw0.euclidianNorm());
		v.add(v0);
		Av.add(Av0);
		
		addQ();
	}
	
	public void addQ()
	{
		final Vector latestQ = latestQNormalized();
		final Vector aq = A.mvMul(latestQ);
		q.add(aq);
	}
	
	public Vector latestQNormalized()
	{
		final Vector lq = q.get(q.size() - 1);
		return lq
			.mul(1. / lq.euclidianNorm());
	}
	
	public Vector latestQ()
	{
		final Vector lq = q.get(q.size() - 1);
		return new DenseVector(lq);
	}
	
	public Vector latestAv()
	{
		return Av.get(v.size() - 1);
	}
	
	public Vector latestV()
	{
		return v.get(v.size() - 1);
	}
	
	public boolean arnoldiStep()
	{
		Vector w = latestQNormalized();
		addQ();
		Vector newAw = latestQ();
		for (int i = 0; i < v.size(); i++)
		{
			final Vector oldV = v.get(i);
			final Vector oldAv = Av.get(i);
			final double projection = newAw.inner(oldAv);
			w = w.sub(oldV.mul(projection));
			newAw = newAw.sub(oldAv.mul(projection));
		}
		final Vector newV = w.mul(1. / newAw.euclidianNorm());
		final Vector newAv = newAw.mul(1. / newAw.euclidianNorm());
		v.add(newV);
		Av.add(newAv);
		//System.out.println("      " + newAw.euclidianNorm());
		return newAw.euclidianNorm() < 1e-13;
	}
}
