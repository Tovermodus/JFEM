package linalg;

import org.junit.Test;

import java.util.HashMap;
import java.util.Map;
import java.util.Random;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

public class BlockDenseMatrixTest
{
	private static BlockDenseMatrix createSmallMatrix()
	{
		final Map<IntCoordinates, DenseMatrix> blocks = new HashMap<>();
		blocks.put(new IntCoordinates(0, 0), new DenseMatrix(4, 4));
		blocks.put(new IntCoordinates(4, 0), new DenseMatrix(2, 4));
		blocks.put(new IntCoordinates(4, 4), new DenseMatrix(2, 2));
		for (int i = 0; i < 16; i++)
			blocks.get(new IntCoordinates(0, 0))
			      .add(i, i / 4, i % 4);
		blocks.get(new IntCoordinates(4, 4))
		      .add(2, 0, 0);
		blocks.get(new IntCoordinates(4, 4))
		      .add(2, 1, 1);
		blocks.get(new IntCoordinates(4, 0))
		      .add(3, 0, 0);
		blocks.get(new IntCoordinates(4, 0))
		      .add(3, 1, 1);
		blocks.get(new IntCoordinates(4, 0))
		      .add(3, 0, 2);
		blocks.get(new IntCoordinates(4, 0))
		      .add(3, 1, 3);
		return new BlockDenseMatrix(blocks);
	}
	
	private static BlockDenseMatrix createMediumMatrix()
	{
		final Map<IntCoordinates, DenseMatrix> blocks = new HashMap<>();
		blocks.put(new IntCoordinates(0, 0), new DenseMatrix(100, 100));
		final BlockDenseMatrix m = new BlockDenseMatrix(blocks);
		for (int i = 0; i < 100 * 100; i++)
			blocks.get(new IntCoordinates(0, 0))
			      .add(2.3 * i, i / 100, i % 100);
		return m;
	}
	
	private static BlockDenseMatrix createLargeMatrix()
	{
		
		final Map<IntCoordinates, DenseMatrix> blocks = new HashMap<>();
		blocks.put(new IntCoordinates(0, 0), new DenseMatrix(150, 150));
		blocks.put(new IntCoordinates(150, 150), new DenseMatrix(150, 150));
		blocks.put(new IntCoordinates(0, 150), new DenseMatrix(150, 150));
		blocks.put(new IntCoordinates(150, 0), new DenseMatrix(150, 150));
		final Random generator = new Random(31415);
		for (final DenseMatrix s : blocks.values())
		{
			for (int i = 0; i < s.getRows(); i++)
				for (int j = 0; j < s.getCols(); j++)
					if (generator.nextDouble() < 0.05)
						s.add(generator.nextDouble() - 0.7, i, j);
		}
		return new BlockDenseMatrix(blocks);
	}
	
	private static BlockDenseMatrix createSmallMatrix2()
	{
		final Map<IntCoordinates, DenseMatrix> blocks = new HashMap<>();
		blocks.put(new IntCoordinates(0, 0), new DenseMatrix(2, 2));
		blocks.put(new IntCoordinates(2, 0), new DenseMatrix(4, 2));
		blocks.put(new IntCoordinates(2, 2), new DenseMatrix(4, 4));
		for (int i = 0; i < 16; i++)
			blocks.get(new IntCoordinates(2, 2))
			      .add(0.01 * i, i / 4, i % 4);
		blocks.get(new IntCoordinates(2, 2))
		      .add(7, 0, 0);
		blocks.get(new IntCoordinates(2, 2))
		      .add(7, 1, 1);
		blocks.get(new IntCoordinates(2, 0))
		      .add(8, 0, 0);
		blocks.get(new IntCoordinates(2, 0))
		      .add(8, 1, 1);
		blocks.get(new IntCoordinates(2, 0))
		      .add(8, 2, 0);
		blocks.get(new IntCoordinates(2, 0))
		      .add(8, 1, 1);
		return new BlockDenseMatrix(blocks);
	}
	
	private static BlockDenseMatrix createMediumMatrix2()
	{
		
		final Map<IntCoordinates, DenseMatrix> blocks = new HashMap<>();
		blocks.put(new IntCoordinates(0, 0), new DenseMatrix(50, 50));
		blocks.put(new IntCoordinates(50, 50), new DenseMatrix(25, 25));
		blocks.put(new IntCoordinates(75, 75), new DenseMatrix(25, 25));
		blocks.put(new IntCoordinates(0, 50), new DenseMatrix(50, 25));
		blocks.put(new IntCoordinates(75, 50), new DenseMatrix(25, 25));
		final Random generator = new Random(3141);
		for (final DenseMatrix s : blocks.values())
		{
			for (int i = 0; i < s.getRows(); i++)
				for (int j = 0; j < s.getCols(); j++)
					s.add(generator.nextDouble(), i, j);
		}
		return new BlockDenseMatrix(blocks);
	}
	
	private static BlockDenseMatrix createLargeMatrix2()
	{
		final Map<IntCoordinates, DenseMatrix> blocks = new HashMap<>();
		blocks.put(new IntCoordinates(0, 0), new DenseMatrix(150, 150));
		blocks.put(new IntCoordinates(150, 150), new DenseMatrix(75, 75));
		blocks.put(new IntCoordinates(225, 225), new DenseMatrix(75, 75));
		blocks.put(new IntCoordinates(0, 150), new DenseMatrix(150, 75));
		blocks.put(new IntCoordinates(225, 150), new DenseMatrix(75, 75));
		final Random generator = new Random(314);
		for (final DenseMatrix s : blocks.values())
		{
			for (int i = 0; i < s.getRows(); i++)
				for (int j = 0; j < s.getCols(); j++)
					if (generator.nextDouble() < 0.05)
						s.add(generator.nextDouble() * 1e7, i, j);
		}
		return new BlockDenseMatrix(blocks);
	}
	
	@Test
	public void testGetters()
	{
		assertTrue(true);
	}
	
	@Test
	public void testFromDense()
	{
		DenseMatrix s = createMediumMatrix().toDense();
		assertEquals(s, new BlockDenseMatrix(s, 1));
		assertEquals(s, new BlockDenseMatrix(s, 2));
		assertEquals(s, new BlockDenseMatrix(s, 3));
		assertEquals(s, new BlockDenseMatrix(s, 8));
		s = createMediumMatrix2().toDense();
		assertEquals(s, new BlockDenseMatrix(s, 1));
		assertEquals(s, new BlockDenseMatrix(s, 2));
		assertEquals(s, new BlockDenseMatrix(s, 3));
		assertEquals(s, new BlockDenseMatrix(s, 8));
		s = createLargeMatrix().toDense();
		assertEquals(s, new BlockDenseMatrix(s, 1));
		assertEquals(s, new BlockDenseMatrix(s, 2));
		assertEquals(s, new BlockDenseMatrix(s, 3));
		assertEquals(s, new BlockDenseMatrix(s, 8));
		assertEquals(s, new BlockDenseMatrix(s, 20));
		assertEquals(s, new BlockDenseMatrix(s, 29));
		s = createLargeMatrix2().toDense();
		assertEquals(s, new BlockDenseMatrix(s, 1));
		assertEquals(s, new BlockDenseMatrix(s, 2));
		assertEquals(s, new BlockDenseMatrix(s, 3));
		assertEquals(s, new BlockDenseMatrix(s, 8));
		assertEquals(s, new BlockDenseMatrix(s, 20));
		assertEquals(s, new BlockDenseMatrix(s, 29));
	}
	
	@Test
	public void testAt()
	{
		final BlockDenseMatrix small = createSmallMatrix();
		for (int j = 0; j < 16; j++)
		{
			assertEquals(j, (int) small.at(j / 4, j % 4));
		}
		assertEquals((int) small.at(4, 0), 3);
		assertEquals((int) small.at(5, 0), 0);
		assertEquals((int) small.at(4, 1), 0);
		assertEquals((int) small.at(5, 1), 3);
		assertEquals((int) small.at(4, 2), 3);
		assertEquals((int) small.at(5, 2), 0);
		assertEquals((int) small.at(4, 3), 0);
		assertEquals((int) small.at(5, 3), 3);
		assertEquals((int) small.at(0, 4), 0);
		assertEquals((int) small.at(0, 5), 0);
		assertEquals((int) small.at(1, 4), 0);
		assertEquals((int) small.at(1, 5), 0);
		assertEquals((int) small.at(2, 4), 0);
		assertEquals((int) small.at(2, 5), 0);
		assertEquals((int) small.at(3, 4), 0);
		assertEquals((int) small.at(3, 5), 0);
		assertEquals((int) small.at(4, 4), 2);
		assertEquals((int) small.at(4, 5), 0);
		assertEquals((int) small.at(5, 4), 0);
		assertEquals((int) small.at(5, 5), 2);
	}
	
	@Test
	public void testAdd()
	{
		final BlockDenseMatrix large = createLargeMatrix();
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
	public void testToDense()
	{
		final DenseMatrix small = new DenseMatrix(6, 6);
		for (int i = 0; i < 16; i++)
			small.add(i, i / 4, i % 4);
		small.add(2, 4, 4);
		small.add(2, 5, 5);
		small.add(3, 4, 0);
		small.add(3, 5, 1);
		small.add(3, 4, 2);
		small.add(3, 5, 3);
		assertEquals(small, createSmallMatrix().toDense());
		assertTrue(small.almostEqual(createSmallMatrix()));
		assertTrue(createSmallMatrix().almostEqual(small));
		assertTrue(createSmallMatrix().toDense()
		                              .almostEqual(createSmallMatrix()));
		assertTrue(createSmallMatrix().almostEqual(createSmallMatrix().toDense()));
		assertTrue(createSmallMatrix().toDense()
		                              .almostEqual(createSmallMatrix().toDense()));
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
		final DenseVector small = new DenseVector(6);
		final DenseVector medium = new DenseVector(100);
		final DenseVector large = new DenseVector(300);
		for (int i = 0; i < 300; i++)
		{
			if (i < 6)
				small.add(i, i);
			if (i < 100)
				medium.add(Math.random(), i);
			large.add(Math.random(), i);
		}
		assertEquals(createSmallMatrix().mvMul(small),
		             createSmallMatrix().toDense()
		                                .mvMul(small));
		assertEquals(createMediumMatrix().mvMul(medium),
		             createMediumMatrix().toDense()
		                                 .mvMul(medium));
		assertEquals(createLargeMatrix().mvMul(large),
		             createLargeMatrix().toDense()
		                                .mvMul(large));
		assertEquals(createSmallMatrix2().mvMul(small),
		             createSmallMatrix2().toDense()
		                                 .mvMul(small));
		assertEquals(createMediumMatrix2().mvMul(medium),
		             createMediumMatrix2().toDense()
		                                  .mvMul(medium));
		assertEquals(createLargeMatrix2().mvMul(large),
		             createLargeMatrix2().toDense()
		                                 .mvMul(large));
	}
	
	@Test
	public void testTvMul()
	{
		final DenseVector small = new DenseVector(6);
		final DenseVector medium = new DenseVector(100);
		final DenseVector large = new DenseVector(300);
		for (int i = 0; i < 300; i++)
		{
			if (i < 6)
				small.add(i, i);
			if (i < 100)
				medium.add(Math.random(), i);
			large.add(Math.random(), i);
		}
		assertEquals(createSmallMatrix().tvMul(small),
		             createSmallMatrix().toDense()
		                                .tvMul(small));
		assertEquals(createMediumMatrix().tvMul(medium),
		             createMediumMatrix().toDense()
		                                 .tvMul(medium));
		assertEquals(createLargeMatrix().tvMul(large),
		             createLargeMatrix().toDense()
		                                .tvMul(large));
		assertEquals(createSmallMatrix2().tvMul(small),
		             createSmallMatrix2().toDense()
		                                 .tvMul(small));
		assertEquals(createMediumMatrix2().tvMul(medium),
		             createMediumMatrix2().toDense()
		                                  .tvMul(medium));
		assertEquals(createLargeMatrix2().tvMul(large),
		             createLargeMatrix2().toDense()
		                                 .tvMul(large));
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
		final BlockDenseMatrix large = createLargeMatrix();
		final DenseMatrix largeSp = large.toDense();
		final BlockDenseMatrix large2 = createLargeMatrix2();
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
