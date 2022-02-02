package linalg;

import java.util.stream.IntStream;

public abstract class ExplicitSchurSolver
	extends AbstractSchurSolver<DenseMatrix>
	implements VectorMultiplyable
{
	final DenseMatrix schurComplement;
	
	public ExplicitSchurSolver(final BlockSparseMatrix blockMatrix)
	{
		super(blockMatrix);
		System.out.println("mmul");
		schurComplement = new DenseMatrix(IntStream.range(0, diagonalInverses.size())
		                                           .parallel()
		                                           .mapToObj(i -> getTopBlock(i)
			                                           .mmMul(diagonalInverses.get(i))
			                                           .mmMul(getLeftBlock(i)))
		                                           .reduce(DenseMatrix::add)
		                                           .orElseThrow());
		System.out.println("-1");
		schurComplement.mulInPlace(-1);
		System.out.println("add");
		schurComplement.addInPlace(getSchurBlock());
	}
	
	public DenseMatrix getSchurComplement()
	{
		return schurComplement;
	}
	
	@Override
	public DenseMatrix schurMvMul()
	{
		return schurComplement;
	}
}
