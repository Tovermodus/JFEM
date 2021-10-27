package linalg;

import org.junit.Test;

import java.util.HashMap;
import java.util.Map;
import java.util.Random;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

public class BlockSparseMatrixTest
{
	private static BlockSparseMatrix createSmallMatrix()
	{
		final Map<IntCoordinates, SparseMatrix> blocks = new HashMap<>();
		blocks.put(new IntCoordinates(0, 0), new SparseMatrix(4, 4));
		blocks.put(new IntCoordinates(4, 0), new SparseMatrix(2, 4));
		blocks.put(new IntCoordinates(4, 4), new SparseMatrix(2, 2));
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
		return new BlockSparseMatrix(blocks);
	}
	
	private static BlockSparseMatrix createMediumMatrix()
	{
		final Map<IntCoordinates, SparseMatrix> blocks = new HashMap<>();
		blocks.put(new IntCoordinates(0, 0), new SparseMatrix(100, 100));
		final BlockSparseMatrix m = new BlockSparseMatrix(blocks);
		for (int i = 0; i < 100 * 100; i++)
			blocks.get(new IntCoordinates(0, 0))
			      .add(2.3 * i, i / 100, i % 100);
		return m;
	}
	
	private static BlockSparseMatrix createLargeMatrix()
	{
		
		final Map<IntCoordinates, SparseMatrix> blocks = new HashMap<>();
		blocks.put(new IntCoordinates(0, 0), new SparseMatrix(500, 500));
		blocks.put(new IntCoordinates(500, 500), new SparseMatrix(500, 500));
		blocks.put(new IntCoordinates(0, 500), new SparseMatrix(500, 500));
		blocks.put(new IntCoordinates(500, 0), new SparseMatrix(500, 500));
		final Random generator = new Random(31415);
		for (final SparseMatrix s : blocks.values())
		{
			for (int i = 0; i < s.getRows(); i++)
				for (int j = 0; j < s.getCols(); j++)
					if (generator.nextDouble() < 0.05)
						s.add(generator.nextDouble(), i, j);
		}
		return new BlockSparseMatrix(blocks);
	}
	
	private static BlockSparseMatrix createSmallMatrix2()
	{
		final Map<IntCoordinates, SparseMatrix> blocks = new HashMap<>();
		blocks.put(new IntCoordinates(0, 0), new SparseMatrix(2, 2));
		blocks.put(new IntCoordinates(2, 0), new SparseMatrix(4, 2));
		blocks.put(new IntCoordinates(2, 2), new SparseMatrix(4, 4));
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
		return new BlockSparseMatrix(blocks);
	}
	
	private static BlockSparseMatrix createMediumMatrix2()
	{
		
		final Map<IntCoordinates, SparseMatrix> blocks = new HashMap<>();
		blocks.put(new IntCoordinates(0, 0), new SparseMatrix(50, 50));
		blocks.put(new IntCoordinates(50, 50), new SparseMatrix(25, 25));
		blocks.put(new IntCoordinates(75, 75), new SparseMatrix(25, 25));
		blocks.put(new IntCoordinates(0, 50), new SparseMatrix(50, 25));
		blocks.put(new IntCoordinates(75, 50), new SparseMatrix(25, 25));
		final Random generator = new Random(3141);
		for (final SparseMatrix s : blocks.values())
		{
			for (int i = 0; i < s.getRows(); i++)
				for (int j = 0; j < s.getCols(); j++)
					s.add(generator.nextDouble(), i, j);
		}
		return new BlockSparseMatrix(blocks);
	}
	
	private static BlockSparseMatrix createLargeMatrix2()
	{
		final Map<IntCoordinates, SparseMatrix> blocks = new HashMap<>();
		blocks.put(new IntCoordinates(0, 0), new SparseMatrix(500, 500));
		blocks.put(new IntCoordinates(500, 500), new SparseMatrix(250, 250));
		blocks.put(new IntCoordinates(750, 750), new SparseMatrix(250, 250));
		blocks.put(new IntCoordinates(0, 500), new SparseMatrix(500, 250));
		blocks.put(new IntCoordinates(750, 500), new SparseMatrix(250, 250));
		final Random generator = new Random(314);
		for (final SparseMatrix s : blocks.values())
		{
			for (int i = 0; i < s.getRows(); i++)
				for (int j = 0; j < s.getCols(); j++)
					if (generator.nextDouble() < 0.05)
						s.add(generator.nextDouble() * 1e7, i, j);
		}
		return new BlockSparseMatrix(blocks);
	}
	
	@Test
	public void testGetters()
	{
		assertTrue(true);
	}
	
	@Test
	public void testAt()
	{
		final BlockSparseMatrix small = createSmallMatrix();
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
		final BlockSparseMatrix large = createLargeMatrix();
		assertEquals(createSmallMatrix().add(createSmallMatrix()),
		             createSmallMatrix().toSparse()
		                                .add(createSmallMatrix().toSparse()));
		assertEquals(createSmallMatrix().add(createSmallMatrix().toSparse()),
		             createSmallMatrix().toSparse()
		                                .add(createSmallMatrix().toSparse()));
		assertEquals(createSmallMatrix().add(createSmallMatrix()),
		             createSmallMatrix().toSparse()
		                                .add(createSmallMatrix()));
		
		assertEquals(createMediumMatrix().add(createMediumMatrix()),
		             createMediumMatrix().toSparse()
		                                 .add(createMediumMatrix().toSparse()));
		assertEquals(createMediumMatrix().add(createMediumMatrix().toSparse()),
		             createMediumMatrix().toSparse()
		                                 .add(createMediumMatrix().toSparse()));
		assertEquals(createMediumMatrix().add(createMediumMatrix()),
		             createMediumMatrix().toSparse()
		                                 .add(createMediumMatrix()));
		
		assertEquals(large.add(large),
		             large.toSparse()
		                  .add(large.toSparse()));
		assertEquals(large.add(large.toSparse()),
		             large.toSparse()
		                  .add(large.toSparse()));
		assertEquals(large.add(large),
		             large.toSparse()
		                  .add(large));
		
		assertEquals(createSmallMatrix2().add(createSmallMatrix()),
		             createSmallMatrix2().toSparse()
		                                 .add(createSmallMatrix().toSparse()));
		assertEquals(createSmallMatrix2().add(createSmallMatrix().toSparse()),
		             createSmallMatrix2().toSparse()
		                                 .add(createSmallMatrix().toSparse()));
		assertEquals(createSmallMatrix2().add(createSmallMatrix()),
		             createSmallMatrix2().toSparse()
		                                 .add(createSmallMatrix()));
		
		assertEquals(createMediumMatrix2().add(createMediumMatrix()),
		             createMediumMatrix2().toSparse()
		                                  .add(createMediumMatrix().toSparse()));
		assertEquals(createMediumMatrix2().add(createMediumMatrix().toSparse()),
		             createMediumMatrix2().toSparse()
		                                  .add(createMediumMatrix().toSparse()));
		assertEquals(createMediumMatrix2().add(createMediumMatrix()),
		             createMediumMatrix2().toSparse()
		                                  .add(createMediumMatrix()));
		
		assertEquals(createLargeMatrix2().add(large),
		             createLargeMatrix2().toSparse()
		                                 .add(large.toSparse()));
		assertEquals(createLargeMatrix2().add(large.toSparse()),
		             createLargeMatrix2().toSparse()
		                                 .add(large.toSparse()));
		assertEquals(createLargeMatrix2().add(large),
		             createLargeMatrix2().toSparse()
		                                 .add(large));
	}
	
	@Test
	public void testMul()
	{
		assertEquals(createSmallMatrix().mul(3)
		                                .toSparse(),
		             createSmallMatrix().toSparse()
		                                .mul(3));
		assertEquals(createSmallMatrix().mul(1), createSmallMatrix());
		assertEquals(createSmallMatrix().add(createSmallMatrix()), createSmallMatrix().mul(2));
		
		assertEquals(createMediumMatrix().mul(3.99999)
		                                 .toSparse(),
		             createMediumMatrix().toSparse()
		                                 .mul(3.99999));
		assertEquals(createMediumMatrix().mul(1), createMediumMatrix());
		assertEquals(createMediumMatrix().add(createMediumMatrix()), createMediumMatrix().mul(2));
		
		assertEquals(createLargeMatrix().mul(1), createLargeMatrix());
		assertEquals(createLargeMatrix().mul(3)
		                                .toSparse(),
		             createLargeMatrix().toSparse()
		                                .mul(3));
		assertEquals(createLargeMatrix().add(createLargeMatrix()), createLargeMatrix().mul(2));
	}
	
	@Test
	public void testToSparse()
	{
		final SparseMatrix small = new SparseMatrix(6, 6);
		for (int i = 0; i < 16; i++)
			small.add(i, i / 4, i % 4);
		small.add(2, 4, 4);
		small.add(2, 5, 5);
		small.add(3, 4, 0);
		small.add(3, 5, 1);
		small.add(3, 4, 2);
		small.add(3, 5, 3);
		assertEquals(small, createSmallMatrix().toSparse());
		assertTrue(small.almostEqual(createSmallMatrix()));
		assertTrue(createSmallMatrix().almostEqual(small));
		assertTrue(createSmallMatrix().toSparse()
		                              .almostEqual(createSmallMatrix()));
		assertTrue(createSmallMatrix().almostEqual(createSmallMatrix().toSparse()));
		assertTrue(createSmallMatrix().toSparse()
		                              .almostEqual(createSmallMatrix().toSparse()));
	}
	
	@Test
	public void testUnfoldDimension()
	{
		assertEquals(createSmallMatrix().unfoldDimension(0),
		             createSmallMatrix().toSparse()
		                                .unfoldDimension(0));
		assertEquals(createSmallMatrix().unfoldDimension(1),
		             createSmallMatrix().toSparse()
		                                .unfoldDimension(1));
		
		assertEquals(createMediumMatrix().unfoldDimension(0),
		             createMediumMatrix().toSparse()
		                                 .unfoldDimension(0));
		assertEquals(createMediumMatrix().unfoldDimension(1),
		             createMediumMatrix().toSparse()
		                                 .unfoldDimension(1));
		
		assertEquals(createLargeMatrix().unfoldDimension(0),
		             createLargeMatrix().toSparse()
		                                .unfoldDimension(0));
		assertEquals(createLargeMatrix().unfoldDimension(1),
		             createLargeMatrix().toSparse()
		                                .unfoldDimension(1));
	}
	
	@Test
	public void testTranspose()
	{
		assertEquals(createSmallMatrix().transpose()
		                                .toSparse(),
		             createSmallMatrix().toSparse()
		                                .transpose());
		assertEquals(createSmallMatrix().transpose()
		                                .transpose(), createSmallMatrix());
		
		assertEquals(createMediumMatrix().transpose()
		                                 .toSparse(),
		             createMediumMatrix().toSparse()
		                                 .transpose());
		assertEquals(createMediumMatrix().transpose()
		                                 .transpose(), createMediumMatrix());
		
		assertEquals(createLargeMatrix().transpose()
		                                .toSparse(),
		             createLargeMatrix().toSparse()
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
		             createSmallMatrix().toSparse()
		                                .mmMul(createSmallMatrix().toSparse()));
		assertEquals(createSmallMatrix().mmMul(createSmallMatrix().toSparse()),
		             createSmallMatrix().toSparse()
		                                .mmMul(createSmallMatrix().toSparse()));
		assertEquals(createSmallMatrix().mmMul(createSmallMatrix()),
		             createSmallMatrix().toSparse()
		                                .mmMul(createSmallMatrix()));
		
		assertEquals(createMediumMatrix().mmMul(createMediumMatrix()),
		             createMediumMatrix().toSparse()
		                                 .mmMul(createMediumMatrix().toSparse()));
		assertEquals(createMediumMatrix().mmMul(createMediumMatrix().toSparse()),
		             createMediumMatrix().toSparse()
		                                 .mmMul(createMediumMatrix().toSparse()));
		assertEquals(createMediumMatrix().mmMul(createMediumMatrix()),
		             createMediumMatrix().toSparse()
		                                 .mmMul(createMediumMatrix()));
		final BlockSparseMatrix large = createLargeMatrix();
		final SparseMatrix largeSp = large.toSparse();
		final BlockSparseMatrix large2 = createLargeMatrix2();
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
		             createSmallMatrix2().toSparse()
		                                 .mmMul(createSmallMatrix().toSparse()));
		assertEquals(createSmallMatrix2().mmMul(createSmallMatrix().toSparse()),
		             createSmallMatrix2().toSparse()
		                                 .mmMul(createSmallMatrix().toSparse()));
		assertEquals(createSmallMatrix2().mmMul(createSmallMatrix()),
		             createSmallMatrix2().toSparse()
		                                 .mmMul(createSmallMatrix()));
		
		assertEquals(createMediumMatrix2().mmMul(createMediumMatrix()),
		             createMediumMatrix2().toSparse()
		                                  .mmMul(createMediumMatrix().toSparse()));
		assertEquals(createMediumMatrix2().mmMul(createMediumMatrix().toSparse()),
		             createMediumMatrix2().toSparse()
		                                  .mmMul(createMediumMatrix().toSparse()));
		assertEquals(createMediumMatrix2().mmMul(createMediumMatrix()),
		             createMediumMatrix2().toSparse()
		                                  .mmMul(createMediumMatrix()));
		
		assertEquals(large2.mmMul(large),
		             large2.toSparse()
		                   .mmMul(largeSp));
		assertEquals(large2.mmMul(largeSp),
		             large2.toSparse()
		                   .mmMul(largeSp));
		assertEquals(large2.mmMul(large),
		             large2.toSparse()
		                   .mmMul(large));
	}
}
