package linalg;

import linalg.DoubleTensor;
import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.LowerTriangPackMatrix;
import no.uib.cipr.matrix.Vector;
import no.uib.cipr.matrix.sparse.LinkedSparseMatrix;
import org.ojalgo.matrix.store.Primitive64Store;
import org.ojalgo.matrix.store.SparseStore;
import org.ujmp.core.DenseMatrix;
import org.ujmp.core.Matrix;
import org.ujmp.core.SparseMatrix;

import javax.swing.*;

interface MatrixFormats
{
	static Matrix getUJMPMatrix(DoubleTensor t)
	{
		Matrix mat;
		if(t.isSparse())
			mat = SparseMatrix.Factory.zeros(t.getRows(),t.getCols());
		else
			mat = DenseMatrix.Factory.zeros(t.getRows(),t.getCols());
		if(!t.isSparse())
			for(int i = 0; i < t.getM(); i++)
				for(int j = 0; j < t.getN(); j++)
					mat.setAsDouble(t.denseValues[i*t.getN()+j], i, j);
		else
			for(int i = 0; i < t.sparseValues.length; i++)
			{
				mat.setAsDouble(mat.getAsDouble(t.sparseYs[i],t.sparseXs[i])+t.sparseValues[i],
					t.sparseYs[i],
					t.sparseXs[i]);
			}
		return mat;
	}
	static SparseStore<Double> getOJALGOMatrix(DoubleTensor t)
	{

		if(t.getDimension() != 2|| !t.isSparse())
			throw new UnsupportedOperationException();
		SparseStore<Double> m = SparseStore.PRIMITIVE64.make(t.getM(),t.getN());
		for(int i = 0; i < t.sparseValues.length; i++)
		{
			m.add(t.sparseYs[i],t.sparseXs[i],t.sparseValues[i]);
		}
		return m;
	}
	static LinkedSparseMatrix getMTJSparseMatrix(DoubleTensor t)
	{
		if(t.getDimension() != 2|| !t.isSparse())
			throw new UnsupportedOperationException();
		LinkedSparseMatrix m= new LinkedSparseMatrix(t.getM(),t.getN());
		for(int i = 0; i < t.sparseValues.length; i++)
		{
			m.add(t.sparseYs[i],t.sparseXs[i],t.sparseValues[i]);
		}
		return m;
	}
	static no.uib.cipr.matrix.DenseMatrix getMTJDenseMatrix(DoubleTensor t)
	{
		if(t.getDimension() != 2)
			throw new UnsupportedOperationException();
		no.uib.cipr.matrix.DenseMatrix m = new no.uib.cipr.matrix.DenseMatrix(t.getM(),t.getN());
		if(t.isSparse())
			for(int i = 0; i < t.sparseValues.length; i++)
			{
				m.add(t.sparseYs[i],t.sparseXs[i],t.sparseValues[i]);
			}
		else
			for(int i = 0; i < t.getM(); i++)
				for(int j = 0; j < t.getN(); j++)
					m.add(i,j,t.at(i,j));
		return m;
	}
	static Matrix getUJMPVector(DoubleTensor t)
	{
		Matrix mat = DenseMatrix.Factory.zeros(t.getRows(),1);
		for(int i = 0; i < t.size(); i++)
			mat.setAsDouble(t.denseValues[i], i, 0);
		return mat;
	}
	static Primitive64Store getOJALGOVector(DoubleTensor t)
	{
		if(t.getDimension() != 1)
			throw new UnsupportedOperationException();
		Primitive64Store m = Primitive64Store.FACTORY.make(t.size(),1);
		for(int i = 0; i < t.size(); i++)
		{
			m.set(i,t.at(i));
		}
		return m;
	}
	static DenseVector getMTJvector(DoubleTensor t)
	{
		if(t.getDimension() != 1)
			throw new UnsupportedOperationException();
		DenseVector m = new DenseVector((int)t.size());
		for(int i = 0; i < t.size(); i++)
		{
			m.set(i,t.at(i));
		}
		return m;
	}
	static DoubleTensor fromUJMPVector(Matrix vector)
	{
		DoubleTensor ret = new DoubleTensor((int)vector.getRowCount());
		for(int i = 0; i < ret.size(); i++)
		{
			ret.denseValues[i] = vector.getAsDouble(i,0);
		}
		return ret;
	}
	static DoubleTensor fromOJALGOVector(Primitive64Store v)
	{
		DoubleTensor m = new DoubleTensor(v.size());
		for(int i = 0; i < v.size(); i++)
		{
			m.set(i,v.get(i));
		}
		return m;
	}
	static DoubleTensor fromMTJvector(Vector v)
	{
		DoubleTensor m = new DoubleTensor(v.size());
		for(int i = 0; i < v.size(); i++)
		{
			m.set(i,v.get(i));
		}
		return m;
	}
	static DoubleTensor fromUJMPMatrix(Matrix m)
	{
		DoubleTensor ret = new DoubleTensor((int)m.getRowCount(),(int)m.getColumnCount(),false);
		for(int i = 0; i < ret.getRows();i++)
			for(int j = 0; j < ret.getCols(); j++)
				ret.set(i,j,m.getAsDouble(i,j));
		return ret;
	}
	static DoubleTensor fromMTJMatrix(no.uib.cipr.matrix.Matrix m)
	{
		DoubleTensor ret = new DoubleTensor((int)m.numRows(),(int)m.numColumns(),false);
		for(int i = 0; i < ret.getRows();i++)
			for(int j = 0; j < ret.getCols(); j++)
				ret.set(i,j,m.get(i,j));
		return ret;
	}
	static LowerTriangPackMatrix getLowerTriangMTJMatrix(DoubleTensor t)
	{
		if(!t.isSparse()||t.getDimension() != 2)
			throw new UnsupportedOperationException();

		LowerTriangPackMatrix m = new LowerTriangPackMatrix(t.getM());
		System.out.println(t.sparseEntries);
		for(int i = 0; i < t.sparseEntries; i++)
		{
			m.add(t.sparseYs[i],t.sparseXs[i],t.sparseValues[i]);
		}
		return m;
	}
	static void printMTJMatrix(no.uib.cipr.matrix.Matrix m)
	{
		System.out.println();
		for(int i = 0; i < m.numRows(); i++)
		{
			System.out.print("[");
			for (int j = 0; j < m.numColumns(); j++)
			{
				if (Math.abs(m.get(i, j)) <= 1e-14)
					System.out.print("<>........  ");
				else
				{
					if (m.get(i, j) >= 0)
						System.out.print("+");
					System.out.print(String.format("%6.3e", m.get(i, j)) + "  ");
				}
			}
			System.out.println("]");
		}
	}

}
