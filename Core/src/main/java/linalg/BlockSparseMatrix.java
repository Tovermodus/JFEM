package linalg;

import basic.MapCollector;
import basic.PerformanceArguments;
import com.google.common.collect.ImmutableMap;
import io.vavr.Tuple2;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.IntStream;
import java.util.stream.Stream;

public class BlockSparseMatrix
	implements Matrix
{
	final ImmutableMap<IntCoordinates, SparseMatrix> blocks;
	final int[] blockStarts;
	final int[] blockEnds;
	
	public BlockSparseMatrix(final Map<IntCoordinates, SparseMatrix> blocks)
	{
		this.blocks = ImmutableMap.copyOf(blocks);
		blockStarts = blocks.keySet()
		                    .stream()
		                    .mapToInt(ic -> ic.get(0))
		                    .distinct()
		                    .sorted()
		                    .toArray();
		blockEnds = new int[blockStarts.length];
		if (blockEnds.length - 1 >= 0) System.arraycopy(blockStarts, 1, blockEnds, 0, blockEnds.length - 1);
		final IntCoordinates lastStart = new IntCoordinates(blockStarts[blockStarts.length - 1],
		                                                    blockStarts[blockStarts.length - 1]);
		blockEnds[blockEnds.length - 1] = lastStart.add(blocks.get(lastStart)
		                                                      .getShape())
		                                           .get(0);
		blocks.forEach(this::checkIfBlockFits);
	}
	
	private void checkIfBlockFits(final IntCoordinates k, final SparseMatrix v)
	{
		int startsX = -1;
		int startsY = -1;
		for (int i = 0; i < blockStarts.length; i++)
		{
			if (k.get(0) == blockStarts[i])
				startsY = i;
			if (k.get(1) == blockStarts[i])
				startsX = i;
		}
		if (startsX == -1 || startsY == -1)
		{
			throw new IllegalArgumentException(
				"Blocks do not fit on grid, especially not the " +
					"block " +
					"starting at " + k);
		}
		if (v.getRows() != blockEnds[startsY] - blockStarts[startsY] || v.getCols() != blockEnds[startsX] - blockStarts[startsX])
			throw new IllegalArgumentException("Block starting at " + k + " has " +
				                                   " wrong size: " +
				                                   v.getRows() + "!=" + (blockEnds[startsY] - blockStarts[startsY]) + " or " + v
				.getCols() +
				                                   "!=" + (blockEnds[startsX] - blockStarts[startsX]));
	}
	
	@Override
	public double at(final int... coordinates)
	{
		if (PerformanceArguments.getInstance().executeChecks)
		{
			if (coordinates.length != 2)
				throw new IllegalArgumentException("Wrong number of coordinates");
			if (coordinates[0] < 0 || coordinates[0] >= getRows())
				throw new IllegalArgumentException("y coordinate out of bounds" + coordinates[0] + " " +
					                                   "with rows: " + getRows());
			if (coordinates[1] < 0 || coordinates[1] >= getCols())
				throw new IllegalArgumentException("x coordinate out of bounds");
		}
		int blockY = 0;
		int blockX = 0;
		for (int i = 0; i < blockStarts.length; i++)
		{
			if (coordinates[0] >= blockStarts[i])
				blockY = i;
			if (coordinates[1] >= blockStarts[i])
				blockX = i;
		}
		final SparseMatrix block = blocks.get(new IntCoordinates(blockStarts[blockY], blockStarts[blockX]));
		if (block == null)
			return 0;
		else
			return block.at(coordinates[0] - blockStarts[blockY],
			                coordinates[1] - blockStarts[blockX]);
	}
	
	@Override
	public Matrix add(final Tensor other)
	{
		return this.toSparse()
		           .add(other);
	}
	
	public SparseMatrix toSparse()
	{
		final SparseMatrix ret = new SparseMatrix(getShape().get(0), getShape().get(1));
		blocks.forEach((key, value) ->
			               value.getCoordinateEntryList()
			                    .forEach((coord, val) ->
				                             ret.add(val, coord.add(key))));
		return ret;
	}
	
	@Override
	public BlockSparseMatrix mul(final double scalar)
	{
		final Map<IntCoordinates, SparseMatrix> ret = new HashMap<>();
		for (final Map.Entry<IntCoordinates, SparseMatrix> entry : blocks.entrySet())
			ret.put(entry.getKey(), entry.getValue()
			                             .mul(scalar));
		return new BlockSparseMatrix(ret);
	}
	
	@Override
	public List<Vector> unfoldDimension(final int dimension)
	{
		return this.toSparse()
		           .unfoldDimension(dimension);
	}
	
	@Override
	public int getSparseEntryCount()
	{
		return blocks.values()
		             .stream()
		             .mapToInt(b -> getSparseEntryCount())
		             .sum();
	}
	
	@Override
	public boolean isSparse()
	{
		return true;
	}
	
	@Override
	public IntCoordinates getShape()
	{
		
		return new IntCoordinates(blockEnds[blockEnds.length - 1], blockEnds[blockEnds.length - 1]);
	}
	
	@Override
	public int getRows()
	{
		return blockEnds[blockEnds.length - 1];
	}
	
	@Override
	public int getCols()
	{
		return blockEnds[blockEnds.length - 1];
	}
	
	@Override
	public boolean almostEqual(final Tensor other)
	{
		if (PerformanceArguments.getInstance().executeChecks)
		{
			if (!getShape().equals(other.getShape()))
				throw new IllegalArgumentException("Tensors are of different size");
		}
		if (!other.isSparse())
			return other.almostEqual(this);
		final ImmutableMap<IntCoordinates, Double> myValues = getCoordinateEntryList();
		final ImmutableMap<IntCoordinates, Double> otherValues = other.getCoordinateEntryList();
		final double absmax = absMaxElement() + other.absMaxElement();
		if (SparseMatrix.coordinateEntryListsEqual(myValues, otherValues, absmax)) return false;
		if (otherValues.size() != myValues.size())
			System.out.println("matrices have different numbers of values");
		return otherValues.size() == myValues.size();
	}
	
	@Override
	public BlockSparseMatrix transpose()
	{
		final Map<IntCoordinates, SparseMatrix> ret = new HashMap<>();
		for (final Map.Entry<IntCoordinates, SparseMatrix> entry : blocks.entrySet())
			ret.put(new IntCoordinates(entry.getKey()
			                                .get(1), entry.getKey()
			                                              .get(0)),
			        entry.getValue()
			             .transpose());
		return new BlockSparseMatrix(ret);
	}
	
	@Override
	public Vector mvMul(final Vector vector)
	{
		final Map<Integer, Vector> subVectors =
			IntStream.range(0, blockStarts.length)
			         .parallel()
			         .mapToObj(i -> new Tuple2<Integer, Vector>(blockStarts[i],
			                                                    vector.slice(blockStarts[i], blockEnds[i])))
			         .collect(new MapCollector<>());
		final DenseVector ret = new DenseVector(vector.getLength());
		final Map<Integer, DenseVector> multipliedSubVectors =
			blocks
				.entrySet()
				.stream()
				.parallel()
				.map(e -> new Tuple2<>(e.getKey()
				                        .get(0), e.getValue()
				                                  .mvMul(subVectors.get(
					                                  e.getKey()
					                                   .get(1)))))
				.collect(new MapCollector<>());
		multipliedSubVectors.forEach((key, value) -> ret.addSmallVectorAt(value, key));
		return ret;
	}
	
	@Override
	public Vector tvMul(final Vector vector)
	{
		
		final Map<Integer, Vector> subVectors =
			IntStream.range(0, blockStarts.length)
			         .parallel()
			         .mapToObj(i -> new Tuple2<Integer, Vector>(blockStarts[i],
			                                                    vector.slice(blockStarts[i], blockEnds[i])))
			         .collect(new MapCollector<>());
		final DenseVector ret = new DenseVector(vector.getLength());
		final Map<Integer, DenseVector> multipliedSubVectors =
			blocks
				.entrySet()
				.stream()
				.parallel()
				.map(e -> new Tuple2<>(e.getKey()
				                        .get(1), e.getValue()
				                                  .tvMul(subVectors.get(
					                                  e.getKey()
					                                   .get(0)))))
				.collect(new MapCollector<>());
		multipliedSubVectors.forEach((key, value) -> ret.addSmallVectorAt(value, key));
		return ret;
	}
	
	@Override
	public Matrix mmMul(final Matrix matrix)
	{
		return toSparse().mmMul(matrix);
	}
	
	@Override
	public String toString()
	{
		return printFormatted();
	}
	
	@Override
	public boolean equals(final Object obj)
	{
		if (!(obj instanceof Matrix))
			return false;
		return almostEqual((Tensor) obj);
	}
	
	@Override
	public ImmutableMap<IntCoordinates, Double> getCoordinateEntryList()
	{
		Stream<Map.Entry<IntCoordinates, SparseMatrix>> str = blocks
			.entrySet()
			.stream();
		if (PerformanceArguments.getInstance().parallelizeThreads)
			str = str.parallel();
		return ImmutableMap.copyOf(str.flatMap(block_entry -> block_entry.getValue()
		                                                                 .getCoordinateEntryList()
		                                                                 .entrySet()
		                                                                 .stream()
		                                                                 .map(coordinate_entry -> new Tuple2<>(
			                                                                 block_entry.getKey()
			                                                                            .add(coordinate_entry
				                                                                                 .getKey()),
			                                                                 coordinate_entry.getValue())))
		                              .collect(new MapCollector<>()));
	}
}
