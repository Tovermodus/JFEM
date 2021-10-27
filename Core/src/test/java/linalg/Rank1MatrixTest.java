package linalg;

import org.junit.Test;

import java.util.Random;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

public class Rank1MatrixTest
{
	private static Rank1Matrix createSmallMatrix()
	{
		return new Rank1Matrix(DenseVector.vectorFromValues(1, 2, 4), DenseVector.vectorFromValues(3, 5, 7));
	}
	
	private static Rank1Matrix createMediumMatrix()
	{
		final DenseVector d1 = new DenseVector(100);
		final DenseVector d2 = new DenseVector(100);
		for (int i = 0; i < 100; i++)
		{
			d1.set(i, i);
			d2.set(i % 17 * 0.29, i);
		}
		return new Rank1Matrix(d1, d2);
	}
	
	private static Rank1Matrix createLargeMatrix()
	{
		final Random generator = new Random(31415);
		final DenseVector d1 = new DenseVector(400);
		final DenseVector d2 = new DenseVector(400);
		for (int i = 0; i < 400; i++)
		{
			d1.set(i, i);
			d2.set(generator.nextDouble(), i);
		}
		return new Rank1Matrix(d1, d2);
	}
	
	private static Rank1Matrix createSmallMatrix2()
	{
		return new Rank1Matrix(DenseVector.vectorFromValues(2, 8, 4), DenseVector.vectorFromValues(-3, 6, 19));
	}
	
	private static Rank1Matrix createMediumMatrix2()
	{
		final Random generator = new Random(3145);
		final DenseVector d1 = new DenseVector(100);
		final DenseVector d2 = new DenseVector(100);
		for (int i = 0; i < 100; i++)
		{
			d1.set(i, i);
			d2.set(generator.nextDouble(), i);
		}
		return new Rank1Matrix(d1, d2);
	}
	
	private static Rank1Matrix createLargeMatrix2()
	{
		final Random generator = new Random(3115);
		final DenseVector d1 = new DenseVector(400);
		final DenseVector d2 = new DenseVector(400);
		for (int i = 0; i < 400; i++)
		{
			d1.set(i, i);
			d2.set(generator.nextDouble(), i);
		}
		return new Rank1Matrix(d1, d2);
	}
	
	@Test
	public void testGetters()
	{
		assertTrue(true);
	}
	
	@Test
	public void testAt()
	{
		final Rank1Matrix small = createSmallMatrix();
	}
	
	@Test
	public void testAdd()
	{
		final Rank1Matrix large = createLargeMatrix();
		assertEquals(createSmallMatrix().add(createSmallMatrix()),
		             createSmallMatrix().toDense()
		                                .add(createSmallMatrix().toDense()));
		assertEquals(createSmallMatrix().add(createSmallMatrix().toDense()),
		             createSmallMatrix().toDense()
		                                .add(createSmallMatrix().toDense()));
		assertEquals(createSmallMatrix().add(createSmallMatrix()),
		             createSmallMatrix().toDense()
		                                .add(createSmallMatrix()));
		
		assertEquals(createMediumMatrix().add(createMediumMatrix()),
		             createMediumMatrix().toDense()
		                                 .add(createMediumMatrix().toDense()));
		assertEquals(createMediumMatrix().add(createMediumMatrix().toDense()),
		             createMediumMatrix().toDense()
		                                 .add(createMediumMatrix().toDense()));
		assertEquals(createMediumMatrix().add(createMediumMatrix()),
		             createMediumMatrix().toDense()
		                                 .add(createMediumMatrix()));
		
		assertEquals(large.add(large),
		             large.toDense()
		                  .add(large.toDense()));
		assertEquals(large.add(large.toDense()),
		             large.toDense()
		                  .add(large.toDense()));
		assertEquals(large.add(large),
		             large.toDense()
		                  .add(large));
		
		assertEquals(createSmallMatrix2().add(createSmallMatrix()),
		             createSmallMatrix2().toDense()
		                                 .add(createSmallMatrix().toDense()));
		assertEquals(createSmallMatrix2().add(createSmallMatrix().toDense()),
		             createSmallMatrix2().toDense()
		                                 .add(createSmallMatrix().toDense()));
		assertEquals(createSmallMatrix2().add(createSmallMatrix()),
		             createSmallMatrix2().toDense()
		                                 .add(createSmallMatrix()));
		
		assertEquals(createMediumMatrix2().add(createMediumMatrix()),
		             createMediumMatrix2().toDense()
		                                  .add(createMediumMatrix().toDense()));
		assertEquals(createMediumMatrix2().add(createMediumMatrix().toDense()),
		             createMediumMatrix2().toDense()
		                                  .add(createMediumMatrix().toDense()));
		assertEquals(createMediumMatrix2().add(createMediumMatrix()),
		             createMediumMatrix2().toDense()
		                                  .add(createMediumMatrix()));
		
		assertEquals(createLargeMatrix2().add(large),
		             createLargeMatrix2().toDense()
		                                 .add(large.toDense()));
		assertEquals(createLargeMatrix2().add(large.toDense()),
		             createLargeMatrix2().toDense()
		                                 .add(large.toDense()));
		assertEquals(createLargeMatrix2().add(large),
		             createLargeMatrix2().toDense()
		                                 .add(large));
	}
	
	@Test
	public void testMul()
	{
		assertEquals(createSmallMatrix().mul(3)
		                                .toDense(),
		             createSmallMatrix().toDense()
		                                .mul(3));
		assertEquals(createSmallMatrix().mul(1), createSmallMatrix());
		assertEquals(createSmallMatrix().add(createSmallMatrix()), createSmallMatrix().mul(2));
		
		assertEquals(createMediumMatrix().mul(3.99999)
		                                 .toDense(),
		             createMediumMatrix().toDense()
		                                 .mul(3.99999));
		assertEquals(createMediumMatrix().mul(1), createMediumMatrix());
		assertEquals(createMediumMatrix().add(createMediumMatrix()), createMediumMatrix().mul(2));
		
		assertEquals(createLargeMatrix().mul(1), createLargeMatrix());
		assertEquals(createLargeMatrix().mul(3)
		                                .toDense(),
		             createLargeMatrix().toDense()
		                                .mul(3));
		assertEquals(createLargeMatrix().add(createLargeMatrix()), createLargeMatrix().mul(2));
	}
	
	@Test
	public void testUnfoldDimension()
	{
		assertEquals(createSmallMatrix().unfoldDimension(0),
		             createSmallMatrix().toDense()
		                                .unfoldDimension(0));
		assertEquals(createSmallMatrix().unfoldDimension(1),
		             createSmallMatrix().toDense()
		                                .unfoldDimension(1));
		
		assertEquals(createMediumMatrix().unfoldDimension(0),
		             createMediumMatrix().toDense()
		                                 .unfoldDimension(0));
		assertEquals(createMediumMatrix().unfoldDimension(1),
		             createMediumMatrix().toDense()
		                                 .unfoldDimension(1));
		
		assertEquals(createLargeMatrix().unfoldDimension(0),
		             createLargeMatrix().toDense()
		                                .unfoldDimension(0));
		assertEquals(createLargeMatrix().unfoldDimension(1),
		             createLargeMatrix().toDense()
		                                .unfoldDimension(1));
	}
	
	@Test
	public void testTranspose()
	{
		assertEquals(createSmallMatrix().transpose()
		                                .toDense(),
		             createSmallMatrix().toDense()
		                                .transpose());
		assertEquals(createSmallMatrix().transpose()
		                                .transpose(), createSmallMatrix());
		
		assertEquals(createMediumMatrix().transpose()
		                                 .toDense(),
		             createMediumMatrix().toDense()
		                                 .transpose());
		assertEquals(createMediumMatrix().transpose()
		                                 .transpose(), createMediumMatrix());
		
		assertEquals(createLargeMatrix().transpose()
		                                .toDense(),
		             createLargeMatrix().toDense()
		                                .transpose());
		assertEquals(createLargeMatrix().transpose()
		                                .transpose(), createLargeMatrix());
	}
	
	@Test
	public void testMvMul()
	{
	}
	
	@Test
	public void testTvMul()
	{
	}
	
	@Test
	public void testMmMul()
	{
		assertEquals(createSmallMatrix().mmMul(createSmallMatrix()),
		             createSmallMatrix().toDense()
		                                .mmMul(createSmallMatrix().toDense()));
		assertEquals(createSmallMatrix().mmMul(createSmallMatrix().toDense()),
		             createSmallMatrix().toDense()
		                                .mmMul(createSmallMatrix().toDense()));
		assertEquals(createSmallMatrix().mmMul(createSmallMatrix()),
		             createSmallMatrix().toDense()
		                                .mmMul(createSmallMatrix()));
		
		assertEquals(createMediumMatrix().mmMul(createMediumMatrix()),
		             createMediumMatrix().toDense()
		                                 .mmMul(createMediumMatrix().toDense()));
		assertEquals(createMediumMatrix().mmMul(createMediumMatrix().toDense()),
		             createMediumMatrix().toDense()
		                                 .mmMul(createMediumMatrix().toDense()));
		assertEquals(createMediumMatrix().mmMul(createMediumMatrix()),
		             createMediumMatrix().toDense()
		                                 .mmMul(createMediumMatrix()));
		final Rank1Matrix large = createLargeMatrix();
		final DenseMatrix largeSp = large.toDense();
		final Rank1Matrix large2 = createLargeMatrix2();
		assertEquals(large.mmMul(large),
		             largeSp
			             .mmMul(largeSp));
		assertEquals(large.mmMul(largeSp),
		             largeSp
			             .mmMul(largeSp));
		assertEquals(large.mmMul(large),
		             largeSp
			             .mmMul(large));
		
		assertEquals(createSmallMatrix2().mmMul(createSmallMatrix()),
		             createSmallMatrix2().toDense()
		                                 .mmMul(createSmallMatrix().toDense()));
		assertEquals(createSmallMatrix2().mmMul(createSmallMatrix().toDense()),
		             createSmallMatrix2().toDense()
		                                 .mmMul(createSmallMatrix().toDense()));
		assertEquals(createSmallMatrix2().mmMul(createSmallMatrix()),
		             createSmallMatrix2().toDense()
		                                 .mmMul(createSmallMatrix()));
		
		assertEquals(createMediumMatrix2().mmMul(createMediumMatrix()),
		             createMediumMatrix2().toDense()
		                                  .mmMul(createMediumMatrix().toDense()));
		assertEquals(createMediumMatrix2().mmMul(createMediumMatrix().toDense()),
		             createMediumMatrix2().toDense()
		                                  .mmMul(createMediumMatrix().toDense()));
		assertEquals(createMediumMatrix2().mmMul(createMediumMatrix()),
		             createMediumMatrix2().toDense()
		                                  .mmMul(createMediumMatrix()));
		
		assertEquals(large2.mmMul(large),
		             large2.toDense()
		                   .mmMul(largeSp));
		assertEquals(large2.mmMul(largeSp),
		             large2.toDense()
		                   .mmMul(largeSp));
		assertEquals(large2.mmMul(large),
		             large2.toDense()
		                   .mmMul(large));
	}
}
