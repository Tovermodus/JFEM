package linalg;

import scala.Function2;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public abstract class SchurSolver
	implements VectorMultiplyable
{
	final BlockSparseMatrix blockMatrix;
	final ArrayList<DenseMatrix> diagonalInverses;
	final ArrayList<SparseMatrix> topMatrices;
	final ArrayList<SparseMatrix> leftMatrices;
	final DenseMatrix schurComplement;
	final int[] blockStarts;
	final int[] blockEnds;
	final int[] blockSizes;
	
	public SchurSolver(final BlockSparseMatrix blockMatrix)
	{
		this.blockMatrix = blockMatrix;
		diagonalInverses = new ArrayList<>();
		topMatrices = new ArrayList<>();
		leftMatrices = new ArrayList<>();
		blockStarts = blockMatrix.getBlockStarts();
		blockEnds = blockMatrix.getBlockEnds();
		blockSizes = blockMatrix.getBlockSizes();
		for (int i = 1; i < blockStarts.length; i++)
		{
			diagonalInverses.add(blockMatrix.getBlocks()
			                                .get(new IntCoordinates(blockStarts[i],
			                                                        blockStarts[i]))
			                                .inverse());
			topMatrices.add(blockMatrix.getBlocks()
			                           .get(new IntCoordinates(0,
			                                                   blockStarts[i])));
			leftMatrices.add(blockMatrix.getBlocks()
			                            .get(new IntCoordinates(blockStarts[i], 0)));
			if (diagonalInverses.get(i - 1) == null)
				throw new IllegalArgumentException("matrix does not contain necessary diagonal block");
			if (topMatrices.get(i - 1) == null)
				throw new IllegalArgumentException("matrix does not contain necessary top block");
			if (leftMatrices.get(i - 1) == null)
				throw new IllegalArgumentException("matrix does not contain necessary left block");
		}
		schurComplement = new DenseMatrix(IntStream.range(0, diagonalInverses.size())
		                                           .parallel()
		                                           .mapToObj(i -> topMatrices.get(i)
		                                                                     .mmMul(diagonalInverses.get(i))
		                                                                     .mmMul(leftMatrices.get(i)))
		                                           .reduce(DenseMatrix::add)
		                                           .orElseThrow());
		schurComplement.mulInPlace(-1);
		schurComplement.addInPlace(blockMatrix.getBlockMatrix(0, 0));
	}
	
	abstract Function2<DenseMatrix, Vector, Vector> solveSchur();
	
	@Override
	public int getVectorSize()
	{
		return blockMatrix.getVectorSize();
	}
	
	@Override
	public int getTVectorSize()
	{
		return blockMatrix.getTVectorSize();
	}
	
	@Override
	public Vector mvMul(final Vector vector)
	{
		final List<DenseVector> subVectors = BlockSparseMatrix.partitionVector(vector,
		                                                                       blockStarts,
		                                                                       blockEnds);
		final DenseVector schurVector = subVectors.get(0);
		subVectors.remove(0);
		IntStream.range(0, diagonalInverses.size())
		         .parallel()
		         .mapToObj(i -> topMatrices.get(i)
		                                   .mvMul(diagonalInverses.get(i)
		                                                          .mvMul(subVectors.get(i))))
		         .forEach(schurVector::subInPlace);
		final List<DenseVector> solvedSubVectors =
			IntStream.range(0, diagonalInverses.size())
			         .parallel()
			         .mapToObj(i -> diagonalInverses.get(i)
			                                        .mvMul(subVectors.get(i)))
			         .collect(
				         Collectors.toList());
		final DenseVector solvedSchurVector = new DenseVector(solveSchur().apply(schurComplement, schurVector));
		final List<DenseVector> newSubVectors
			= IntStream.range(0, diagonalInverses.size())
			           .parallel()
			           .mapToObj(i -> solvedSubVectors.get(i)
			                                          .sub(diagonalInverses.get(i)
			                                                               .mvMul(leftMatrices.get(i)
			                                                                                  .mvMul(solvedSchurVector))))
			           .collect(Collectors.toList());
		final DenseVector ret = new DenseVector(getVectorSize());
		ret.addSmallVectorAt(solvedSchurVector, 0);
		IntStream.range(0, diagonalInverses.size())
		         .forEach(i -> ret.addSmallVectorAt(newSubVectors.get(i),
		                                            blockStarts[i + 1]));
		return ret;
	}
	
	@Override
	public Vector tvMul(final Vector vector)
	{
		throw new UnsupportedOperationException("not implemented yet");
	}
}
