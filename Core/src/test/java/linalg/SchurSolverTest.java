package linalg;

import org.junit.Test;

import java.util.HashMap;
import java.util.Map;
import java.util.Random;

import static org.junit.Assert.assertEquals;

public class SchurSolverTest
{
	private static BlockSparseMatrix createSmallIdentity3()
	{
		final Map<IntCoordinates, SparseMatrix> blocks = new HashMap<>();
		blocks.put(new IntCoordinates(0, 0), SparseMatrix.identity(4));
		blocks.put(new IntCoordinates(4, 4), SparseMatrix.identity(2));
		blocks.put(new IntCoordinates(6, 6), SparseMatrix.identity(14));
		blocks.put(new IntCoordinates(0, 4), new SparseMatrix(4, 2));
		blocks.put(new IntCoordinates(4, 0), new SparseMatrix(2, 4));
		blocks.put(new IntCoordinates(0, 6), new SparseMatrix(4, 14));
		blocks.put(new IntCoordinates(6, 0), new SparseMatrix(14, 4));
		return new BlockSparseMatrix(blocks, 20, 20);
	}
	
	private static BlockSparseMatrix createSmallMatrix3()
	{
		final Random generator = new Random(3145);
		final Map<IntCoordinates, SparseMatrix> blocks = new HashMap<>();
		blocks.put(new IntCoordinates(0, 0), SparseMatrix.identity(4));
		blocks.put(new IntCoordinates(4, 4), SparseMatrix.identity(2));
		blocks.put(new IntCoordinates(6, 6), SparseMatrix.identity(14));
		blocks.put(new IntCoordinates(0, 4), new SparseMatrix(4, 2));
		blocks.put(new IntCoordinates(4, 0), new SparseMatrix(2, 4));
		blocks.put(new IntCoordinates(0, 6), new SparseMatrix(4, 14));
		blocks.put(new IntCoordinates(6, 0), new SparseMatrix(14, 4));
		blocks.forEach((k, mat) ->
		               {
			               for (final IntCoordinates c : mat.getShape()
			                                                .range())
			               {
				               mat.add(generator.nextDouble() * 0.02, c);
			               }
		               });
		return new BlockSparseMatrix(blocks, 20, 20);
	}
	
	private static BlockSparseMatrix createLargeMatrix3()
	{
		final Random generator = new Random(38145);
		final Map<IntCoordinates, SparseMatrix> blocks = new HashMap<>();
		blocks.put(new IntCoordinates(0, 0), SparseMatrix.identity(100));
		blocks.put(new IntCoordinates(100, 100), SparseMatrix.identity(50));
		blocks.put(new IntCoordinates(150, 150), SparseMatrix.identity(150));
		blocks.put(new IntCoordinates(0, 100), new SparseMatrix(100, 50));
		blocks.put(new IntCoordinates(100, 0), new SparseMatrix(50, 100));
		blocks.put(new IntCoordinates(0, 150), new SparseMatrix(100, 150));
		blocks.put(new IntCoordinates(150, 0), new SparseMatrix(150, 100));
		blocks.forEach((k, mat) ->
		               {
			               for (final IntCoordinates c : mat.getShape()
			                                                .range())
			               {
				               mat.add(generator.nextDouble() * 0.002, c);
			               }
		               });
		return new BlockSparseMatrix(blocks, 300, 300);
	}
	
	private static BlockSparseMatrix createSmallIdentity2()
	{
		final Map<IntCoordinates, SparseMatrix> blocks = new HashMap<>();
		blocks.put(new IntCoordinates(0, 0), SparseMatrix.identity(4));
		blocks.put(new IntCoordinates(4, 4), SparseMatrix.identity(2));
		blocks.put(new IntCoordinates(0, 4), new SparseMatrix(4, 2));
		blocks.put(new IntCoordinates(4, 0), new SparseMatrix(2, 4));
		return new BlockSparseMatrix(blocks, 6, 6);
	}
	
	private static BlockSparseMatrix createSmallMatrix2()
	{
		final Random generator = new Random(31415);
		final Map<IntCoordinates, SparseMatrix> blocks = new HashMap<>();
		blocks.put(new IntCoordinates(0, 0),
		           SparseMatrix.identity(2)
		                       .mul(1000));
		blocks.put(new IntCoordinates(2, 2),
		           SparseMatrix.identity(1)
		                       .mul(1000));
		blocks.put(new IntCoordinates(0, 2), new SparseMatrix(2, 1));
		blocks.put(new IntCoordinates(2, 0), new SparseMatrix(1, 2));
		blocks.forEach((k, mat) ->
		               {
			               for (final IntCoordinates c : mat.getShape()
			                                                .range())
			               {
				               mat.add(generator.nextInt(10) - 1, c);
			               }
		               });
		return new BlockSparseMatrix(blocks, 3, 3);
	}
	
	private static DenseVector getRhs(final Matrix m)
	{
		final Random generator = new Random(3415);
		final DenseVector d = new DenseVector(m.getVectorSize());
		for (final IntCoordinates c : d.getShape()
		                               .range())
		{
			d.add(generator.nextInt(10), c);
		}
		return d;
	}
	
	@Test
	public void testIdentity()
	{
		final BlockSparseMatrix bsm = createSmallIdentity2();
		final DenseVector b = getRhs(bsm);
		final SchurSolver s = new DirectSchur(bsm);
		assertEquals(s.mvMul(b), b);
		assertEquals(s.mvMul(b),
		             bsm.toSparse()
		                .solve(b));
		final BlockSparseMatrix bsm3 = createSmallIdentity3();
		final DenseVector b3 = getRhs(bsm3);
		final SchurSolver s3 = new DirectSchur(bsm3);
		assertEquals(s3.mvMul(b3), b3);
		assertEquals(s3.mvMul(b3),
		             bsm3.toSparse()
		                 .solve(b3));
	}
	
	@Test
	public void testMatrix()
	{
		final BlockSparseMatrix bsm = createSmallMatrix2();
		final DenseVector b = getRhs(bsm);
		final SchurSolver s = new DirectSchur(bsm);
		assertEquals(s.mvMul(b),
		             bsm.toSparse()
		                .solve(b));
		final BlockSparseMatrix bsm3 = createSmallMatrix3();
		final DenseVector b3 = getRhs(bsm3);
		final SchurSolver s3 = new DirectSchur(bsm3);
		assertEquals(s3.mvMul(b3),
		             bsm3.toSparse()
		                 .solve(b3));
		final BlockSparseMatrix bsml3 = createLargeMatrix3();
		final DenseVector bl3 = getRhs(bsml3);
		final SchurSolver sl3 = new DirectSchur(bsml3);
		assertEquals(sl3.mvMul(bl3),
		             bsml3.toSparse()
		                  .solve(bl3));
	}
}
