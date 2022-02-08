package linalg;

import java.util.stream.IntStream;

public abstract class ExplicitSchurSolver
	extends AbstractSchurSolver<SparseMatrix>
	implements VectorMultiplyable
{
	SparseMatrix schurComplement;
	
	public ExplicitSchurSolver(final BlockSparseMatrix blockMatrix)
	{
		super(blockMatrix);
		System.out.println("mmul");
		schurComplement = new SparseMatrix(IntStream.range(0, diagonalInverses.size())
		                                            .parallel()
		                                            .mapToObj(i -> getTopBlock(i)
			                                            .mmMul(new SparseMatrix(diagonalInverses.get(i)))
			                                            .mmMul(getLeftBlock(i)))
		                                            .reduce(SparseMatrix::add)
		                                            .orElseThrow());
		System.out.println("-1");
		schurComplement.mulInPlace(-1);
		System.out.println("add");
		schurComplement.addInPlace(getSchurBlock());
	}
	
	@Override
	public void resetOffDiagonals(final BlockSparseMatrix systemMatrix)
	{
		super.resetOffDiagonals(systemMatrix);
		System.out.println("mmul");
		schurComplement = new SparseMatrix(IntStream.range(0, diagonalInverses.size())
		                                            .parallel()
		                                            .mapToObj(i -> getTopBlock(i)
			                                            .mmMul(new SparseMatrix(diagonalInverses.get(i)))
			                                            .mmMul(getLeftBlock(i)))
		                                            .reduce(SparseMatrix::add)
		                                            .orElseThrow());
		System.out.println("-1");
		schurComplement.mulInPlace(-1);
		System.out.println("add");
		schurComplement.addInPlace(getSchurBlock());
	}
	
	public SparseMatrix getSchurComplement()
	{
		return schurComplement;
	}
	
	@Override
	public SparseMatrix schurMvMul()
	{
		return schurComplement;
	}
}
