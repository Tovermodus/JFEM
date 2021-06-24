package linalg;

import org.junit.jupiter.api.Test;
import java.util.*;

import static org.junit.jupiter.api.Assertions.*;
public class IntCoordinatesTest
{
	@Test
	public void testGetCoordinates()
	{
		IntCoordinates a1 = new IntCoordinates(1,2,3,4);
		IntCoordinates a2 = new IntCoordinates(0);
		IntCoordinates a3 = new IntCoordinates(1,2,3,4,5,6,7,8,9,0);
		assertArrayEquals(a1.asArray(), new int[]{1, 2, 3, 4});
		assertArrayEquals(a2.asArray(), new int[]{0});
		assertArrayEquals(a3.asArray(), new int[]{1,2,3,4,5,6,7,8,9,0});
	}
	@Test
	public void testGetDimension() {
		IntCoordinates a1 = new IntCoordinates(1,2,3,4);
		IntCoordinates a2 = new IntCoordinates(0);
		IntCoordinates a3 = new IntCoordinates(1,2,3,4,5,6,7,8,9,0);
		assertEquals(a1.getDimension(), 4);
		assertEquals(a2.getDimension(), 1);
		assertEquals(a3.getDimension(), 10);
	}
	@Test
	public void testGetSize() {
		IntCoordinates a1 = new IntCoordinates(1,2,3,4);
		IntCoordinates a2 = new IntCoordinates(2,3,4);
		assertSame(a1.size(), a2.size());
		assertEquals(2*3*4, a1.size());
	}
	@Test
	public void testHashCode() {
		IntCoordinates a1 = new IntCoordinates(1,2,3,4);
		IntCoordinates a2 = new IntCoordinates(0);
		IntCoordinates a3 = new IntCoordinates(1,2,3,4);
		assertEquals(a1.hashCode(), a3.hashCode());
		assertNotEquals(a1.hashCode(), a2.hashCode());
	}
	@Test
	public void testEquals() {
		IntCoordinates a1 = new IntCoordinates(1,2,3,4);
		IntCoordinates a2 = new IntCoordinates(0);
		IntCoordinates a3 = new IntCoordinates(1,2,3,4);
		assertEquals(a1, a3);
		assertNotEquals(a1, a2);
	}
	@Test
	public void testCompare() {
		IntCoordinates a1 = new IntCoordinates(1,2,3,4);
		IntCoordinates a2 = new IntCoordinates(1,2,3);
		IntCoordinates a3 = new IntCoordinates(1,2,3,5);
		IntCoordinates a4 = new IntCoordinates(1,2,3,2);
		IntCoordinates a5 = new IntCoordinates(1,2,6,4);
		IntCoordinates a6 = new IntCoordinates(1,2,1,40);
		IntCoordinates a7 = new IntCoordinates(100,2,3,4);
		IntCoordinates a8 = new IntCoordinates(0,2,3,40);
		IntCoordinates a9 = new IntCoordinates(1,2,3,4);
		assertThrows(IllegalArgumentException.class, () -> a1.compareTo(a2));
		assertEquals(a1.compareTo(a3), -1);
		assertEquals(a3.compareTo(a1), 1);
		assertEquals(a1.compareTo(a4), 1);
		assertEquals(a4.compareTo(a1), -1);
		assertEquals(a1.compareTo(a5), -1);
		assertEquals(a5.compareTo(a1), 1);
		assertEquals(a1.compareTo(a6), 1);
		assertEquals(a6.compareTo(a1), -1);
		assertEquals(a1.compareTo(a7), -1);
		assertEquals(a7.compareTo(a1), 1);
		assertEquals(a1.compareTo(a8), 1);
		assertEquals(a8.compareTo(a1), -1);
		assertEquals(a1.compareTo(a9), 0);
		assertEquals(a9.compareTo(a1), 0);
	}
	@Test
	public void testSet() {
		Set<IntCoordinates> hash = new HashSet<>(10);
		hash.add(new IntCoordinates(1,2,3,4));
		hash.add(new IntCoordinates(0,0,0,0));
		hash.add(new IntCoordinates(0,0,0,0));
		hash.add(new IntCoordinates(7,8,9,0));
		hash.add(new IntCoordinates(7,8,9,0));
		hash.add(new IntCoordinates(1,2,3,4));
		hash.add(new IntCoordinates(7,8,9,0));
		hash.add(new IntCoordinates(1,2,3,4));
		assertEquals(hash.size(), 3);
		Set<IntCoordinates> tree = new TreeSet<>();
		tree.add(new IntCoordinates(1,2,3,4));
		tree.add(new IntCoordinates(0,0,0,0));
		tree.add(new IntCoordinates(0,0,0,0));
		tree.add(new IntCoordinates(7,8,9,0));
		tree.add(new IntCoordinates(7,8,9,0));
		tree.add(new IntCoordinates(1,2,3,4));
		tree.add(new IntCoordinates(7,8,9,0));
		tree.add(new IntCoordinates(1,2,3,4));
		assertEquals(tree.size(), 3);
	}
	@Test
	public void testRangeValid() {
		assertThrows(IllegalArgumentException.class, () -> new IntCoordinates.Range(new int[]{5},new int[]{2}));
		assertThrows(IllegalArgumentException.class, () -> new IntCoordinates.Range(new int[]{5},new int[]{5}));
		assertThrows(IllegalArgumentException.class, () -> new IntCoordinates.Range(new int[]{0,1,2,5},
			new int[]{4,4,4,2}));
		assertThrows(IllegalArgumentException.class, () -> new IntCoordinates.Range(new int[]{0,1,2,1},
			new int[]{0,4,4,2}));
		assertThrows(IllegalArgumentException.class, () -> new IntCoordinates.Range(new int[]{5},new int[]{5,6}));
	}
	@Test
	public void testRange() {
		Set<IntCoordinates> tree = new TreeSet<>();
		IntCoordinates lowerBounds = new IntCoordinates(0,0,1);
		IntCoordinates upperBounds = new IntCoordinates(2,3,5);
		IntCoordinates coordinates = new IntCoordinates(0,0,0);
		for(IntCoordinates c:new IntCoordinates.Range(lowerBounds,upperBounds))
		{
			tree.add(new IntCoordinates(c));
			assertNotEquals(coordinates, c);
			assertEquals(-1, coordinates.compareTo(c));
			coordinates = new IntCoordinates(c);
		}
		assertEquals(2 * 3 * 4, tree.size());
		tree = new TreeSet<>();
		upperBounds = new IntCoordinates(2,3,4);
		coordinates = new IntCoordinates(0,0,0);
		for(IntCoordinates c:new IntCoordinates.Range(upperBounds))
		{
			tree.add(c);
			assertTrue(0 >= coordinates.compareTo(c));
			coordinates = c;
		}
		assertEquals(2 * 3 * 4, tree.size());
	}
	@Test
	public void testRangeStream() {
		final Set<IntCoordinates> tree = new TreeSet<>();
		IntCoordinates lowerBounds = new IntCoordinates(0,0,1);
		IntCoordinates upperBounds = new IntCoordinates(1,3,5);
		new IntCoordinates.Range(lowerBounds,upperBounds).stream().forEach(tree::add);
		assertEquals(1 * 3 * 4, tree.size());
	}
}
