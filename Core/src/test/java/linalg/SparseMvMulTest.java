package linalg;

import com.google.common.base.Stopwatch;
import org.junit.Test;

import java.util.Random;

import static org.junit.Assert.assertEquals;

public class SparseMvMulTest
{
	public static SparseMatrix createTinySparseMatrix()
	{
		final Random generator = new Random(14239078);
		final SparseMatrix ret = new SparseMatrix(10, 10);
		for (final IntCoordinates c : ret.getShape()
		                                 .range())
		{
			if (generator.nextDouble() < 0.1)
				ret.add((int) (10 * generator.nextDouble() - 5), c);
		}
		return ret;
	}
	
	public static SparseMatrix createSmallSparseMatrix()
	{
		final Random generator = new Random(891765234);
		final SparseMatrix ret = new SparseMatrix(500, 700);
		for (final IntCoordinates c : ret.getShape()
		                                 .range())
		{
			if (generator.nextDouble() < 0.01)
				ret.add(generator.nextDouble() - 0.4, c);
		}
		return ret;
	}
	
	public static SparseMatrix createLargeSparseMatrix()
	{
		
		final Random generator = new Random(1627894);
		final SparseMatrix ret = new SparseMatrix(5000, 7000);
		for (final IntCoordinates c : ret.getShape()
		                                 .range())
		{
			if (generator.nextDouble() < 0.05)
				ret.add(generator.nextDouble() - 0.4, c);
		}
		return ret;
	}
	
	public static DenseVector createSmallVector()
	{
		final Random generator = new Random(124876);
		final DenseVector ret = new DenseVector(700);
		for (final IntCoordinates c : ret.getShape()
		                                 .range())
		{
			ret.add(generator.nextDouble() - 0.4, c);
		}
		return ret;
	}
	
	public static DenseVector createTinyVector()
	{
		final Random generator = new Random(124350987);
		final DenseVector ret = new DenseVector(10);
		for (final IntCoordinates c : ret.getShape()
		                                 .range())
		{
			ret.add((int) (10 * generator.nextDouble() - 5), c);
		}
		return ret;
	}
	
	public static DenseVector createLargeVector()
	{
		final Random generator = new Random(18429760);
		
		final DenseVector ret = new DenseVector(7000);
		for (final IntCoordinates c : ret.getShape()
		                                 .range())
		{
			ret.add(generator.nextDouble() - 0.4, c);
		}
		return ret;
	}
	
	public static DenseVector createSmallTVector()
	{
		final Random generator = new Random(81726534);
		
		final DenseVector ret = new DenseVector(500);
		for (final IntCoordinates c : ret.getShape()
		                                 .range())
		{
			ret.add(generator.nextDouble() - 0.4, c);
		}
		return ret;
	}
	
	public static DenseVector createLargeTVector()
	{
		
		final Random generator = new Random(323781645);
		final DenseVector ret = new DenseVector(5000);
		for (final IntCoordinates c : ret.getShape()
		                                 .range())
		{
			ret.add(generator.nextDouble() - 0.4, c);
		}
		return ret;
	}
	
	@Test
	public void getVectorSize()
	{
		assertEquals(new SparseMvMul(createSmallSparseMatrix()).getVectorSize(), 700);
		assertEquals(new SparseMvMul(createLargeSparseMatrix()).getVectorSize(), 7000);
	}
	
	@Test
	public void getTVectorSize()
	{
		assertEquals(new SparseMvMul(createSmallSparseMatrix()).getTVectorSize(), 500);
		assertEquals(new SparseMvMul(createLargeSparseMatrix()).getTVectorSize(), 5000);
	}
	
	@Test
	public void mvMul()
	{
		final SparseMatrix large = createLargeSparseMatrix();
		final SparseMvMul largemv = new SparseMvMul(createLargeSparseMatrix());
		Stopwatch s = Stopwatch.createStarted();
		for (int i = 0; i < 1000; i++)
			large.mvMul(createLargeVector());
		System.out.println(s.elapsed());
		s = Stopwatch.createStarted();
		for (int i = 0; i < 1000; i++)
			largemv.mvMul(createLargeVector());
		System.out.println(s.elapsed());
//		assertEquals(createSmallSparseMatrix().mvMul(createSmallVector()),
//		             new SparseMvMul(createSmallSparseMatrix()).mvMul(createSmallVector()));
//		assertEquals(createTinySparseMatrix().mvMul(createTinyVector()),
//		             new SparseMvMul(createTinySparseMatrix()).mvMul(createTinyVector()));
//		assertEquals(createLargeSparseMatrix().mvMul(createLargeVector()),
//		             new SparseMvMul(createLargeSparseMatrix()).mvMul(createLargeVector()));
	}
	
	@Test
	public void tvMul()
	{
		assertEquals(createTinySparseMatrix().tvMul(createTinyVector()),
		             new SparseMvMul(createTinySparseMatrix()).tvMul(createTinyVector()));
		assertEquals(createSmallSparseMatrix().tvMul(createSmallTVector()),
		             new SparseMvMul(createSmallSparseMatrix()).tvMul(createSmallTVector()));
		assertEquals(createLargeSparseMatrix().tvMul(createLargeTVector()),
		             new SparseMvMul(createLargeSparseMatrix()).tvMul(createLargeTVector()));
	}
}
