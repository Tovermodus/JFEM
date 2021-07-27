package linalg;

import java.util.ArrayList;
import java.util.stream.IntStream;

class GMRES
{
	int restarts = 0;
	public linalg.Vector solve(VectorMultiplyable A,
	                           Vector b,
	                           double tol, IterativeSolver.Interruptor interruptor)
	{
		
		MutableVector x = new linalg.DenseVector(b.getLength());
		x.set(1,0);
		return solve(A,b,x,tol, interruptor);
	}
	public linalg.Vector solve(VectorMultiplyable A,
	                           Vector b, Vector x,
	                           double tol, IterativeSolver.Interruptor interruptor)
	{
		ArrayList<Vector> v = new ArrayList<>();
		int n = b.getLength();
		MutableVector c = new linalg.DenseVector(n);
		MutableVector s = new linalg.DenseVector(n);
		MutableVector gamma = new linalg.DenseVector(n+1);
		Vector r = b.sub(A.mvMul(x));
		SparseMatrix h = new SparseMatrix(n, n);
		if (r.euclidianNorm() <= tol)
			return x;
		gamma.set(r.euclidianNorm(), 0);
		v.add(r.mul(1./gamma.at(0)));
		int j;
		for(j = 0; j < n && r.euclidianNorm() > tol && interruptor.isRunning(); j++)
		{
			Vector q = A.mvMul(v.get(j));
			DenseVector newHValues =
				DenseVector.vectorFromValues(v.stream().parallel().mapToDouble(vec->vec.inner(q)).toArray());
			h.addSmallMatrixAt(newHValues.asMatrix(), 0, j);
			Vector newHV =
				(IntStream.range(0, j+1)).parallel().mapToObj(i->v.get(i).mul(newHValues.at(i))).reduce(new DenseVector(n), Vector::add);
			Vector w = q.sub(newHV);
			h.add(w.euclidianNorm(), j+1,j);
			for(int i = 0; i < j; i++)
			{
				double c_i = c.at(i);
				double s_i = s.at(i);
				double h_ij = h.at(i,j);
				double h_ipj = h.at(i+1,j);
				h.set(h_ij*c_i + h_ipj*s_i, i, j);
				h.set(-h_ij*s_i + h_ipj*c_i, i+1, j);
			}
			double beta = Math.sqrt(Math.pow(h.at(j, j),2) + Math.pow(h.at(j+1,j),2));
			s.set(h.at(j+1,j)/beta, j);
			c.set(h.at(j,j)/beta, j);
			h.set(beta, j, j);
			gamma.set(-s.at(j)*gamma.at(j),j+1);
			gamma.set(c.at(j)*gamma.at(j),j);
			//System.out.println(Math.abs(gamma.at(j+1)));
			if(Math.abs(gamma.at(j+1)) < tol || j > 40)
				break;
			v.add(w.mul(1./h.at(j+1,j)));
		}
		if(j == n || !interruptor.isRunning())
			j--;
		DenseVector alpha = new DenseVector(n);
		for (int i = j; i >=0 ; i--)
		{
			int finalI = i;
			double hAlpha = IntStream.range(i+1,j+1).mapToDouble(k->h.at(finalI,k)*alpha.at(k)).sum();
			alpha.set(1./h.at(i,i)*(gamma.at(i) - hAlpha), i);
		}
		Vector alphaV =
			IntStream.range(0,j+1).parallel().mapToObj(i->v.get(i).mul(alpha.at(i))).reduce(new DenseVector(n),
				Vector::add);
		if(j > 40 && restarts < 200)
		{
			restarts++;
			return new GMRES().solve(A, b, x.add(alphaV), tol, interruptor);
		}
		return x.add(alphaV);
		
	}
	public<T extends VectorMultiplyable> linalg.Vector solve(T preconditioner, VectorMultiplyable A,
	                           Vector b,
	                           double tol, IterativeSolver.Interruptor interruptor)
	{
		
		MutableVector x = new linalg.DenseVector(b.getLength());
		x.set(1,0);
		return solve(preconditioner,A,b,x,tol, interruptor);
	}
	public<T extends VectorMultiplyable> linalg.Vector solve(T preconditioner, VectorMultiplyable A,
	                           Vector b, Vector x,
	                           double tol, IterativeSolver.Interruptor interruptor)
	{
		ArrayList<Vector> v = new ArrayList<>();
		int n = b.getLength();
		MutableVector c = new linalg.DenseVector(n);
		MutableVector s = new linalg.DenseVector(n);
		MutableVector gamma = new linalg.DenseVector(n+1);
		Vector r = preconditioner.mvMul(b.sub(A.mvMul(x)));
		SparseMatrix h = new SparseMatrix(n, n);
		if (r.euclidianNorm() <= tol)
			return x;
		gamma.set(r.euclidianNorm(), 0);
		v.add(r.mul(1./gamma.at(0)));
		int j;
		for(j = 0; j < n && r.euclidianNorm() > tol && interruptor.isRunning(); j++)
		{
			Vector q = preconditioner.mvMul(A.mvMul(v.get(j)));
			DenseVector newHValues =
				DenseVector.vectorFromValues(v.stream().parallel().mapToDouble(vec->vec.inner(q)).toArray());
			h.addSmallMatrixAt(newHValues.asMatrix(), 0, j);
			Vector newHV =
				(IntStream.range(0, j+1)).parallel().mapToObj(i->v.get(i).mul(newHValues.at(i))).reduce(new DenseVector(n), Vector::add);
			Vector w = q.sub(newHV);
			h.add(w.euclidianNorm(), j+1,j);
			for(int i = 0; i < j; i++)
			{
				double c_i = c.at(i);
				double s_i = s.at(i);
				double h_ij = h.at(i,j);
				double h_ipj = h.at(i+1,j);
				h.set(h_ij*c_i + h_ipj*s_i, i, j);
				h.set(-h_ij*s_i + h_ipj*c_i, i+1, j);
			}
			double beta = Math.sqrt(Math.pow(h.at(j, j),2) + Math.pow(h.at(j+1,j),2));
			s.set(h.at(j+1,j)/beta, j);
			c.set(h.at(j,j)/beta, j);
			h.set(beta, j, j);
			gamma.set(-s.at(j)*gamma.at(j),j+1);
			gamma.set(c.at(j)*gamma.at(j),j);
			//System.out.println(Math.abs(gamma.at(j+1)));
			if(Math.abs(gamma.at(j+1)) < tol || j > 40)
				break;
			v.add(w.mul(1./h.at(j+1,j)));
		}
		if(j == n || !interruptor.isRunning())
			j--;
		DenseVector alpha = new DenseVector(n);
		for (int i = j; i >=0 ; i--)
		{
			int finalI = i;
			double hAlpha = IntStream.range(i+1,j+1).mapToDouble(k->h.at(finalI,k)*alpha.at(k)).sum();
			alpha.set(1./h.at(i,i)*(gamma.at(i) - hAlpha), i);
		}
		Vector alphaV =
			IntStream.range(0,j+1).parallel().mapToObj(i->v.get(i).mul(alpha.at(i))).reduce(new DenseVector(n),
				Vector::add);
		if(j > 40 && restarts < 200)
		{
			restarts++;
			return solve(A, preconditioner, b, x.add(alphaV), tol, interruptor);
		}
		return x.add(alphaV);
		
	}
}