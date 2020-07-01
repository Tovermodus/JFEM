package linalg;

import no.uib.cipr.matrix.DenseLU;
import no.uib.cipr.matrix.DenseMatrix;
import no.uib.cipr.matrix.DenseVector;
import no.uib.cipr.matrix.sparse.IterativeSolverNotConvergedException;
import org.jetbrains.annotations.NotNull;
import org.ujmp.core.Matrix;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.stream.IntStream;

public class DoubleTensor implements VectorMultiplication
{
	double [] denseValues;
	double [] sparseValues;
	int [] sparseXs;
	int [] sparseYs;
	int sparseEntries;
	private ArrayList<SparseEntry> sparseEntryArrayList;
	private int cols;
	private int rows;
	private boolean isSparse;
	private int dimension;
	private DenseLU lu;
	public DoubleTensor(int size, double [] denseValues)
	{
		this.denseValues = denseValues;
		this.rows = size;
		this.cols = 1;
		this.isSparse = false;
		this.dimension = 1;
	}
	public DoubleTensor(long size)
	{
		dimension = 1;
		rows = (int)size;
		cols = 1;
		this.isSparse=false;
		denseValues = new double[(int)size];
	}
	public DoubleTensor(int sizeY, int sizeX, boolean isSparse)
	{
		dimension = 2;
		this.isSparse = isSparse;
		if(!isSparse)
			denseValues = new double[sizeX*sizeY];
		this.cols = sizeX;
		this.rows = sizeY;
		if(isSparse)
		{
			sparseValues = new double[1];
			sparseYs = new int[1];
			sparseXs = new int[1];
			this.sparseEntries = 1;
		}
	}
	public DoubleTensor(@NotNull DoubleTensor t)
	{
		cols = t.cols;
		rows = t.rows;
		if(t.isSparse)
		{
			sparseXs = t.sparseXs.clone();
			sparseYs = t.sparseYs.clone();
			sparseValues = t.sparseValues.clone();
		}
		else
			denseValues = t.denseValues.clone();
		isSparse = t.isSparse;
		dimension = t.dimension;
		sparseEntries = t.sparseEntries;
		sparseEntryArrayList = t.sparseEntryArrayList;
	}
	
	public int getDimension()
	{
		return dimension;
	}
	
	public int getCols()
	{
		return cols;
	}
	
	public int getRows()
	{
		return rows;
	}
	
	public boolean isSparse()
	{
		return isSparse;
	}
	
	public void fillSparseEntryList()
	{
		sparseEntryArrayList = new ArrayList<>(sparseEntries);
		for(int i = 0; i < sparseEntries; i++)
		{
			sparseEntryArrayList.add(new SparseEntry(sparseYs[i],sparseXs[i],sparseValues[i]));
		}
		sparseEntryArrayList.sort((entry1,entry2)->(entry1.y - entry2.y)*cols+entry1.x-entry2.x);
	}

	public long size()
	{
		return rows* cols;
	}
	public int getM()
	{
		if(dimension != 2)
			throw new UnsupportedOperationException();
		return rows;
	}
	public int getN()
	{
		if(dimension != 2)
			throw new UnsupportedOperationException();
		return cols;
	}
	@Override
	public boolean equals(Object T)
	{
		if(!(T instanceof DoubleTensor))
			return false;
		if(((DoubleTensor) T).dimension != dimension)
			return false;
		if(((DoubleTensor) T).getRows() != getRows())
			return false;
		if(((DoubleTensor) T).getCols() != getCols())
			return false;
		if(((DoubleTensor) T).dimension == 2)
			for(int i = 0;  i < ((DoubleTensor) T).getRows(); i++)
				for(int j = 0; j < ((DoubleTensor) T).getCols(); j++)
					if(((DoubleTensor) T).at(i,j) != at(i,j))
					{
						System.out.println(((DoubleTensor) T).at(i,j) - at(i,j));
						return false;
					}
		if(((DoubleTensor) T).dimension == 1)
			for(int i = 0;  i < ((DoubleTensor) T).getRows(); i++)
				if(((DoubleTensor) T).at(i) != at(i))
					return false;
		return true;
	}
	public double absmax()
	{
		double ret = 0;
		if(isSparse)
			for(int i = 0; i < sparseEntries; i++)
				if(Math.abs(sparseValues[i]) >ret)
					ret = Math.abs(sparseValues[i]);
		if(!isSparse)
			for(int i = 0; i < denseValues.length; i++)
				if(Math.abs(denseValues[i]) > ret )
						ret = Math.abs(denseValues[i]);
		return ret;
	}
	public boolean almostequals(DoubleTensor T, double tol)
	{
		DoubleTensor n = T.sub(this);
		if(n.dimension == 2)
			for(int i = 0; i < getRows(); i++)
				for(int j = 0; j < getCols(); j++)
					if(Math.abs(n.at(i,j)) > absmax()*tol)
					{
						System.out.println(Math.abs(n.at(i,j))+ " > " + absmax()*tol);
						return false;
					}
		if(n.dimension == 1)
			for(int i = 0; i < getRows(); i++)
				if(Math.abs(n.at(i)) > (Math.abs(this.at(i))+Math.abs(this.at(i)))*tol)
				{
					System.out.println(Math.abs(n.at(i))+ " > " + (Math.abs(this.at(i))+Math.abs(this.at(i)))*tol);
					return false;
				}
		return true;
	}
	public double at(int vectorIndex)
	{
		if(dimension != 1)
			throw new UnsupportedOperationException();
		return denseValues[vectorIndex];
	}
	public double at(int matrixIndex1, int matrixIndex2)
	{
		if(dimension != 2)
			throw new UnsupportedOperationException();
		if(!isSparse)
			return denseValues[matrixIndex1*getN()+matrixIndex2];
		else
		{
			double ret = 0;
			for(int i = 0; i < sparseEntries; i++)
			{
				if(sparseYs[i] == matrixIndex1)
					if(sparseXs[i] == matrixIndex2)
						ret += sparseValues[i];
			}
			return ret;
		}
	}
	public void set(int vectorIndex, double value)
	{
		if(dimension != 1)
			throw new UnsupportedOperationException();
		denseValues[vectorIndex] = value;
	}
	private void resizeSparse()
	{
		double [] sparseVals = new double[sparseValues.length*2];
		int [] sparseX  = new int[sparseValues.length*2];
		int [] sparseY =  new int[sparseValues.length*2];
		for(int i = 0; i < sparseEntries; i++)
		{
			sparseVals[i] = sparseValues[i];
			sparseX[i] = sparseXs[i];
			sparseY[i] = sparseYs[i];
		}
		sparseValues = sparseVals;
		sparseYs = sparseY;
		sparseXs = sparseX;
	}

	synchronized public  void set(int matrixIndex1,int matrixIndex2, double value)
	{
		if(dimension != 2)
			throw new UnsupportedOperationException();
		if(!isSparse)
			denseValues[matrixIndex1*getN()+matrixIndex2] = value;
		else
		{
			if(sparseEntries >= sparseValues.length-2)
				resizeSparse();
			sparseValues[sparseEntries] = value - at(matrixIndex1,matrixIndex2);
			sparseYs[sparseEntries] = matrixIndex1;
			sparseXs[sparseEntries] = matrixIndex2;
			sparseEntries++;
		}

	}
	synchronized public  void add(int matrixIndex1,int matrixIndex2, double value)
	{
		if(dimension != 2)
			throw new UnsupportedOperationException();
		if(isSparse)
		{
			if(sparseEntries >= sparseValues.length-2)
				resizeSparse();
			sparseValues[sparseEntries] = value;
			sparseYs[sparseEntries] = matrixIndex1;
			sparseXs[sparseEntries] = matrixIndex2;
			sparseEntries++;
		}
		else
			denseValues[matrixIndex1*getN()+matrixIndex2]+= value;
	}
	public void add(int vectorIndex, double value)
	{
		if(dimension != 1)
			throw new UnsupportedOperationException();
		denseValues[vectorIndex]+= value;
	}
	public void print_formatted()
	{
		System.out.println();
		if(dimension == 2)
			for(int i = 0; i < rows; i++)
			{
				System.out.print("[");
				for(int j = 0; j < cols; j++)
				{
					if (Math.abs(at(i, j)) <= 1e-14)
						System.out.print("<>........  ");
					else
					{
						if(at(i,j)>=0)
							System.out.print("+");
						System.out.print(String.format("%6.3e", at(i, j)) + "  ");
					}
				}
				System.out.println("]");
			}
		if(dimension == 1)
		{
			System.out.print("[");
			for (int j = 0; j < size(); j++)
			{
				if (Math.abs(at(j)) <= 1e-14)
					System.out.print("<>........  ");
				else
				{
					if(at(j)>=0)
						System.out.print("+");
					System.out.print(String.format("%6.3e", at(j)) + "  ");
				}
			}
			System.out.println("]");
		}
		System.out.println();
	}
	public void print_formatted(int xS,int yS,int xE, int yE)
	{
		System.out.println();
		if(dimension == 2)
			for(int i = yS; i < yE; i++)
			{
				System.out.print("[");
				for(int j = xS; j < xE; j++)
				{
					if (Math.abs(at(i, j)) <= 1e-14)
						System.out.print("<>........  ");
					else
					{
						if(at(i,j)>=0)
							System.out.print("+");
						System.out.print(String.format("%6.3e", at(i, j)) + "  ");
					}
				}
				System.out.println("]");
			}
		if(dimension == 1)
		{
			System.out.print("[");
			for (int j = 0; j < size(); j++)
			{
				if (Math.abs(at(j)) <= 1e-14)
					System.out.print("<>........  ");
				else
				{
					if(at(j)>=0)
						System.out.print("+");
					System.out.print(String.format("%6.3e", at(j)) + "  ");
				}
			}
			System.out.println("]");
		}
		System.out.println();
	}

	public void print_formatted(String s)
	{
		System.out.println("   "+s);
		if(dimension == 2)
			for(int i = 0; i < getM(); i++)
			{
				System.out.print("[");
				for(int j = 0; j < getN(); j++)
				{
					if (Math.abs(at(i, j)) <= 1e-14)
						System.out.print("<>........  ");
					else
					{
						if(at(i,j)>=0)
							System.out.print("+");
						System.out.print(String.format("%6.3e", at(i, j)) + "  ");
					}
				}
				System.out.println("]");
			}
		if(dimension == 1)
		{
			System.out.print("[");
			for (int j = 0; j < size(); j++)
			{
				if (Math.abs(at(j)) <= 1e-14)
					System.out.print("<>........  ");
				else
				{
					if(at(j)>=0)
						System.out.print("+");
					System.out.print(String.format("%6.3e", at(j)) + "  ");
				}
			}
			System.out.println("]");
		}
		System.out.println();
	}
	public void print_log()
	{
		if(dimension != 2)
			return;
		for(int i = 0; i < getM(); i++)
		{
			for(int j = 0; j < getN(); j++)
			{
				if (Math.abs(at(i, j)) <= 1e-14)
					System.out.print("   \t");
				else
					System.out.print((int) Math.log10(Math.abs(at(j, i))) + "\t");
			}
			System.out.println();
		}
	}
	public static DoubleTensor vectorFromValues(double... values)
	{

		DoubleTensor ret = new DoubleTensor(values.length);
		for(int i = 0; i < values.length; i++)
			ret.set(i,values[i]);
		return ret;
	}
	public static DoubleTensor squareMatrixFromValues(double... values)
	{
		if(Math.pow((int)Math.sqrt(values.length),2) != values.length)
			throw new UnsupportedOperationException("nope");
		DoubleTensor ret = new DoubleTensor((int)Math.sqrt(values.length),(int)Math.sqrt(values.length),false);
		ret.denseValues = values;
		return ret;
	}
	public static DoubleTensor identity(int n)
	{
		DoubleTensor ret = new DoubleTensor(n,n,true);
		for(int  i = 0; i < n; i++)
			ret.add(i,i,1);
		return ret;
	}
	public double x()
	{
		if(dimension == 1)
			return at(0);
		else throw new UnsupportedOperationException("nope");
	}
	public double y()
	{
		if(dimension == 1)
			return at(1);
		else throw new UnsupportedOperationException("nope");
	}
	public double inner(@NotNull DoubleTensor other)
	{
		double ret = 0;
		if(this.size() == other.size())
		{

			return IntStream.range(0, (int)this.size())
				.parallel()
				.mapToDouble( id -> denseValues[id] * other.denseValues[id])
				.reduce(0, Double::sum);
		}
		throw new UnsupportedOperationException();
	}
	public DoubleTensor add(@NotNull DoubleTensor other)
	{
		if(this.dimension == other.dimension && this.size() == other.size())
		{
			DoubleTensor ret = null;
			if(!this.isSparse && !other.isSparse)
			{
				ret = new DoubleTensor(this);
				for(int i = 0; i < denseValues.length; i++)
					ret.denseValues[i] = this.denseValues[i]+other.denseValues[i];
			}
			if(!this.isSparse && other.isSparse)
			{
				ret = new DoubleTensor(this);
				for(int i = 0; i < other.sparseEntries; i++)
				{
					ret.denseValues[other.sparseYs[i]*getN() + other.sparseXs[i]] += other.sparseValues[i];
				}
			}
			if(this.isSparse && !other.isSparse)
			{
				ret = new DoubleTensor(other);
				ret.denseValues = other.denseValues.clone();
				for(int i = 0; i < this.sparseEntries; i++)
				{
					ret.denseValues[this.sparseYs[i]*getN() + this.sparseXs[i]] += this.sparseValues[i];
				}
			}
			if(this.isSparse && other.isSparse)
			{
				ret = new DoubleTensor(this.getM(),this.getN(), true);
				for(int i = 0; i < this.sparseEntries; i++)
				{
					ret.add(this.sparseYs[i],this.sparseXs[i],this.sparseValues[i]);
				}
				for(int i = 0; i < other.sparseEntries; i++)
				{
					ret.add(other.sparseYs[i],other.sparseXs[i],other.sparseValues[i]);
				}
			}
			return ret;
		}
		throw new UnsupportedOperationException();
	}
	public DoubleTensor mul(double scalar)
	{
		DoubleTensor ret = new DoubleTensor(this);
		if(!isSparse)
		{
			ret.denseValues = Arrays.stream(denseValues).parallel().map(val->val*scalar).toArray();
		}
		else
		{
			ret.denseValues = Arrays.stream(sparseValues).parallel().map(val->val*scalar).toArray();
		}
		ret.dimension = dimension;
		ret.isSparse = isSparse;
		return ret;
	}
	public DoubleTensor mvmul(@NotNull DoubleTensor vector)
	{
		if(vector.dimension != 1 || this.dimension != 2)
			throw new UnsupportedOperationException();
		DoubleTensor ret = new DoubleTensor(getM());
		ret.isSparse = false;
		ret.dimension = 1;
		if(!isSparse)
		{
			for(int i = 0; i < getM(); i++)
			{
				for(int j = 0; j < getN(); j++)
					ret.denseValues[i] += denseValues[i*getN()+j]*vector.denseValues[j];
			}
		}
		else
		{
			for (int i = 0; i < sparseEntries; i++)
			{
				ret.denseValues[sparseYs[i]] += sparseValues[i] * vector.denseValues[sparseXs[i]];
			}

		}
		return ret;
	}

	public DoubleTensor tvmul(@NotNull DoubleTensor vector)
	{
		if(vector.dimension != 1 || this.dimension != 2)
			throw new UnsupportedOperationException();
		DoubleTensor ret = new DoubleTensor(getN());
		ret.isSparse = false;
		ret.dimension = 1;
		if(!isSparse)
		{
			for(int i = 0; i < getN(); i++)
			{
				for(int j = 0; j < getM(); j++)
					ret.denseValues[i] += denseValues[j*getN()+i]*vector.denseValues[j];
			}
		}
		else
		{
			for (int i = 0; i < sparseEntries; i++)
			{
				ret.denseValues[sparseXs[i]] += sparseValues[i] * vector.denseValues[sparseYs[i]];
			}

		}
		return ret;
	}

	public DoubleTensor mmmul(@NotNull DoubleTensor other)
	{
		if(other.dimension != 2 || this.dimension != 2)
			throw new UnsupportedOperationException();
		DoubleTensor ret = new DoubleTensor(this.getM(),other.getN(),this.isSparse && other.isSparse);
		if(!this.isSparse && !other.isSparse)
		{
			for(int i = 0; i < getM(); i++)
				for(int j = 0; j < other.getN(); j++)
					for(int k = 0; k < getN(); k++)
						ret.denseValues[i*other.getN()+j] +=
							denseValues[i*getN()+k]*other.denseValues[k*other.getN()+j];
		}
		if(!this.isSparse && other.isSparse)
		{
			for(int i = 0; i < other.sparseEntries; i++)
				for(int j = 0; j < getM(); j++)
				{
					ret.denseValues[j * ret.getN() + other.sparseXs[i]] +=
						denseValues[j * getN() + other.sparseYs[i]] * other.sparseValues[i];
				}
		}
		if(this.isSparse && !other.isSparse)
		{
			for(int i = 0; i < this.sparseEntries; i++)
				for(int j = 0; j < other.getN(); j++)
					ret.denseValues[sparseYs[i]*ret.getN() + j] +=
						sparseValues[i] * other.denseValues[sparseXs[i]*other.getN()+j];
		}
		if(this.isSparse && other.isSparse)
		{
			//for(int i = 0; i < this.sparseEntries; i++)
			IntStream.range(0,this.sparseEntries).parallel().forEach(i->
			{
				for (int j = 0; j < other.sparseEntries; j++)
				{
					if (this.sparseXs[i] == other.sparseYs[j])
						ret.add(this.sparseYs[i], other.sparseXs[j],
							this.sparseValues[i] * other.sparseValues[j]);
				}
			});
		}
		return ret;
	}
	public DoubleTensor transpose()
	{
		if(this.dimension != 2)
			throw new UnsupportedOperationException();
		DoubleTensor ret = new DoubleTensor(this);
		ret.dimension = 2;
		if(isSparse)
		{
			ret.sparseXs = sparseYs;
			ret.sparseYs = sparseXs;
			ret.sparseValues = sparseValues;
			ret.sparseEntries = sparseEntries;
		}
		else
		{
			for(int i = 0; i < getM(); i++)
				for(int j = 0; j < getN(); j++)
				{
					ret.denseValues[j*getM()+i] = denseValues[i*getN()+j];
				}
		}
		ret.rows = cols;
		ret.cols = rows;
		ret.isSparse = isSparse;
		return ret;
	}
	public DoubleTensor getLowerTriangleMatrix()
	{
		assert(rows == cols);
		DoubleTensor ret = new DoubleTensor(getRows(), getCols(), isSparse());
		if (isSparse)
		{
			for (int i = 0; i < sparseEntries; i++)
				if (sparseYs[i] < sparseXs[i])
					ret.add(sparseYs[i], sparseXs[i], sparseValues[i]);
		}
		else
			for (int i = 0; i < rows; i++)
				for(int j = 0; j < i; j++)
					ret.add(i,j,at(i,j));
		return ret;
	}
	public DoubleTensor getUpperTriangleMatrix()
	{
		assert(rows == cols);
		DoubleTensor ret = new DoubleTensor(getRows(), getCols(), isSparse());
		if (isSparse)
		{
			for (int i = 0; i < sparseEntries; i++)
				if (sparseYs[i] > sparseXs[i])
					ret.add(sparseYs[i], sparseXs[i], sparseValues[i]);
		}
		else
			for (int i = 0; i < rows; i++)
				for(int j = i+1; j < cols; j++)
					ret.add(i,j,at(i,j));
		return ret;
	}
	public DoubleTensor getDiagonalMatrix()
	{
		assert(rows == cols);
		DoubleTensor ret = new DoubleTensor(getRows(), getCols(), isSparse());
		if (isSparse)
		{
			for (int i = 0; i < sparseEntries; i++)
				if (sparseYs[i] == sparseXs[i])
					ret.add(sparseYs[i], sparseXs[i], sparseValues[i]);
		}
		else
			for (int i = 0; i < rows; i++)
				ret.add(i,i,at(i,i));
		return ret;
	}
	public double vectorNorm()
	{
		if(this.dimension != 1)
			throw new UnsupportedOperationException();
		return Math.sqrt(this.inner(this));
	}
	public DoubleTensor sub(@NotNull DoubleTensor other)
	{
		if(this.dimension != other.dimension && this.size() != other.size())
		{
			throw new UnsupportedOperationException();
		}
		return this.add(other.mul(-1.));
	}
	public DoubleTensor solveSymm(@NotNull DoubleTensor rightHandSide)
	{
		if(rightHandSide.dimension!=1)
			throw new UnsupportedOperationException();
		DoubleTensor ret;
		ret =
			MatrixFormats.fromUJMPVector(MatrixFormats.getUJMPMatrix(this).solveSymm(MatrixFormats.getUJMPVector(rightHandSide)));
		return ret;
	}
	public DoubleTensor inverse()
	{
		Matrix m = MatrixFormats.getUJMPMatrix(this);
		return MatrixFormats.fromUJMPMatrix(m.inv());
	}
	public void generateLU()
	{
		lu = DenseLU.factorize(MatrixFormats.getMTJDenseMatrix(this));
	}
	public DoubleTensor solve(DoubleTensor rightHandSide)
	{
		if(rightHandSide.dimension!=1)
			throw new UnsupportedOperationException();
		if(lu == null)
			generateLU();
		DenseMatrix rhs = new DenseMatrix(MatrixFormats.getMTJvector(rightHandSide));
		DenseMatrix sol = lu.solve(rhs);
		DenseVector sv = new DenseVector(sol.getData());
		return MatrixFormats.fromMTJvector(sv);
	}
	public DoubleTensor solveCG(DoubleTensor rightHandSide, double residualTolerance)
	{
		//System.out.println(this.mat.isSparse());

		if(getN() != getN() || dimension != 2 || rightHandSide.dimension != 1)
			throw new UnsupportedOperationException();
		DoubleTensor z;
		DoubleTensor newResiduum;
		DoubleTensor iterate = new DoubleTensor(getM());
		DoubleTensor residuum = rightHandSide.sub(this.mvmul(iterate));
		DoubleTensor defect = new DoubleTensor(residuum);
		double alpha;
		double beta;
		for(int iter = 0; iter < getM() && residuum.vectorNorm() > residualTolerance; iter++)
		{
			z = this.mvmul(defect);
			alpha = residuum.inner(residuum)/defect.inner(z);
			iterate = iterate.add(defect.mul(alpha));
			newResiduum = residuum.sub(z.mul(alpha));
			beta = newResiduum.inner(newResiduum)/residuum.inner(residuum);
			defect = newResiduum.add(defect.mul(beta));
			residuum = newResiduum;
			if(iter%10 == 0)
				System.out.println("iter "+iter+" " +residuum.vectorNorm());
		}
		return iterate;
	}
	public<T extends VectorMultiplication> DoubleTensor solvePGMRES(T preconditioner, DoubleTensor rightHandSide,
	                                                                double tol)
	{
		DenseVector v = new DenseVector(getM());
		DoubleTensor x = null;
		GMRES gm = new GMRES(v);
		try
		{
			x =  gm.solve(preconditioner, this, rightHandSide, tol);
		} catch (IterativeSolverNotConvergedException e)
		{
			e.printStackTrace();
		}
		return x;
	}
	public DoubleTensor solveGMRES(DoubleTensor rightHandSide, double tol)
	{
		DenseVector v = new DenseVector(getM());
		DoubleTensor x = null;
		GMRES gm = new GMRES(v);
		try
		{
			x =  gm.solve(this, rightHandSide, tol);
		} catch (IterativeSolverNotConvergedException e)
		{
			e.printStackTrace();
		}
		return x;
	}
//	public DoubleTensor solveLowerTriangular(DoubleTensor rightHandSide)
//	{
//		if(this.dimension != 2 || rightHandSide.dimension != 1 || getN()!= getM() || !isSparse)
//			throw new UnsupportedOperationException();
//		DoubleTensor ret = new DoubleTensor(rightHandSide.size());
//		fillSparseEntryList();
//		for(SparseEntry entry:sparseEntryArrayList)
//		{
//
//		}
//
//	}
	public DoubleTensor solveBiCGStab(DoubleTensor rightHandSide, double residualTolerance)
	{
		if(getN() != getN() || dimension != 2 || rightHandSide.dimension != 1)
			throw new UnsupportedOperationException();
		DoubleTensor v;
		DoubleTensor s;
		DoubleTensor t;
		DoubleTensor iterate = new DoubleTensor(getM());
		DoubleTensor startResiduum = rightHandSide.sub(this.mvmul(iterate));
		DoubleTensor residuum = new DoubleTensor(startResiduum);
		DoubleTensor p = new DoubleTensor(residuum);
		double alpha;
		double omega;
		double rho = residuum.inner(residuum);
		double rhoLast;
		double beta;
		for(int iter = 0; iter < getM() && residuum.vectorNorm() > residualTolerance; iter++)
		{
			v = this.mvmul(p);
			alpha = rho/v.inner(startResiduum);
			s = residuum.sub(v.mul(alpha));
			t = this.mvmul(s);
			omega = t.inner(s)/t.inner(t);
			iterate = iterate.add(p.mul(alpha)).add(s.mul(omega));
			residuum = s.sub(t.mul(omega));
			rhoLast = rho;
			rho = residuum.inner(startResiduum);
			beta = alpha/omega*rho/rhoLast;
			p = residuum.add(p.mul(beta)).sub(v.mul(omega*beta));
			if(iter%3 == 0)
				System.out.println("iter "+iter+" "+residuum.vectorNorm());
		}
		return iterate;


	}
	public<T extends VectorMultiplication> DoubleTensor solvePBiCGStab(T preconditioner, DoubleTensor rightHandSide,
	                                                                             double residualTolerance)
	{
		if(getN() != getN() || dimension != 2 || rightHandSide.dimension != 1)
			throw new UnsupportedOperationException();
		DoubleTensor v;
		DoubleTensor vP;
		DoubleTensor s;
		DoubleTensor sP;
		DoubleTensor t;
		DoubleTensor tP;
		DoubleTensor iterate = new DoubleTensor(getM());
		DoubleTensor startResiduum = rightHandSide.sub(this.mvmul(iterate));
		DoubleTensor startResiduumP = preconditioner.mvmul(startResiduum);
		DoubleTensor residuum = new DoubleTensor(startResiduum);
		DoubleTensor residuumP = preconditioner.mvmul(residuum);
		DoubleTensor pP = preconditioner.mvmul(residuumP);
		double alpha;
		double omega;
		double rho = residuumP.inner(residuumP);
		double rhoLast;
		double beta;
		for(int iter = 0; iter < getM() && residuum.vectorNorm() > residualTolerance; iter++)
		{
			v = this.mvmul(pP);
			vP = preconditioner.mvmul(v);
			alpha = rho/vP.inner(startResiduumP);
			s = residuum.sub(v.mul(alpha));
			sP = preconditioner.mvmul(s);
			t = this.mvmul(sP);
			tP = preconditioner.mvmul(t);
			omega = tP.inner(sP)/tP.inner(tP);
			iterate = iterate.add(pP.mul(alpha)).add(sP.mul(omega));
			residuum = rightHandSide.sub(this.mvmul(iterate));
			residuumP = sP.sub(tP.mul(omega));
			rhoLast = rho;
			rho = residuumP.inner(startResiduumP);
			beta = alpha/omega*rho/rhoLast;
			pP = residuumP.add(pP.mul(beta)).sub(vP.mul(omega*beta));
			if(iter%3 == 0)
				System.out.println("iter "+iter+" "+residuum.vectorNorm());
		}
		return iterate;


	}
}

class SparseEntry
{
	public int x;
	public int y;
	public double value;

	public SparseEntry(int sparseY, int sparseX, double sparseValue)
	{
		x = sparseX;
		y = sparseY;
		value = sparseValue;
	}
}
