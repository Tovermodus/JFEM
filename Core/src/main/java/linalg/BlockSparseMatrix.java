package linalg;

import basic.DoubleCompare;
import basic.MapCollector;
import basic.PerformanceArguments;
import com.google.common.collect.ImmutableMap;
import io.vavr.Tuple2;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import java.util.stream.Stream;

public class BlockSparseMatrix
	implements Matrix
{
	public ImmutableMap<IntCoordinates, SparseMatrix> getBlocks()
	{
		return blocks;
	}
	
	public SparseMatrix getBlockMatrix(final IntCoordinates blockNumber)
	{
		return blocks.get(blockNumber.concatenate(blockStarts));
	}
	
	public SparseMatrix getBlockMatrix(final int... blockNumber)
	{
		return getBlockMatrix(new IntCoordinates(blockNumber));
	}
	
	public int[] getBlockStarts()
	{
		return blockStarts.clone();
	}
	
	public int[] getBlockEnds()
	{
		return blockEnds.clone();
	}
	
	public BlockSparseMatrix subMatrix(final IntCoordinates[] blockNumbers)
	{
		final Map<IntCoordinates, SparseMatrix> newBlocks = new HashMap<>();
		for (int i = 0; i < blockStarts.length; i++)
			newBlocks.put(new IntCoordinates(blockStarts[i], blockStarts[i]),
			              new SparseMatrix(blockEnds[i] - blockStarts[i], blockEnds[i] - blockStarts[i]));
		for (final IntCoordinates blockNumber : blockNumbers)
		{
			newBlocks.put(blockNumber.concatenate(blockStarts), getBlockMatrix(blockNumber));
		}
		return new BlockSparseMatrix(newBlocks);
	}
	
	final private ImmutableMap<IntCoordinates, SparseMatrix> blocks;
	final private int[] blockStarts;
	final private int[] blockEnds;
	
	public BlockSparseMatrix(final Map<IntCoordinates, SparseMatrix> blocks)
	{
		this.blocks = ImmutableMap.copyOf(blocks);
		blockStarts = blocks.keySet()
		                    .stream()
		                    .mapToInt(c -> c.get(0))
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
	
	public BlockSparseMatrix(final SparseMatrix s, final int nBlocksPerDir)
	{
		this(s, generateEquiDist(nBlocksPerDir, s.getRows()));
	}
	
	private static int[] generateEquiDist(final int nBlocksPerDir, final int rows)
	{
		final int[] blockst = new int[nBlocksPerDir];
		for (int i = 0; i < nBlocksPerDir; i++)
			blockst[i] = (int) (1.0 * i * rows / nBlocksPerDir);
		return blockst;
	}
	
	public BlockSparseMatrix(final SparseMatrix s, final int[] blockStarts)
	{
		if (PerformanceArguments.getInstance().executeChecks)
		{
			if (s.getRows() != s.getCols())
				throw new IllegalArgumentException("Must be square");
		}
		this.blockStarts = blockStarts;
		blockEnds = new int[blockStarts.length];
		if (blockEnds.length - 1 >= 0) System.arraycopy(blockStarts, 1, blockEnds, 0, blockEnds.length - 1);
		blockEnds[blockEnds.length - 1] = s.getCols();
		final HashMap<IntCoordinates, SparseMatrix> mutableBlocks = new HashMap<>();
		for (final Map.Entry<IntCoordinates, Double> entry : s.getCoordinateEntryList()
		                                                      .entrySet())
		{
			final IntCoordinates block = getBlock(entry.getKey()
			                                           .asArray());
			final IntCoordinates blockStart = block.concatenate(blockStarts);
			if (!mutableBlocks.containsKey(blockStart))
				mutableBlocks.put(blockStart,
				                  new SparseMatrix(block.concatenate(blockEnds)
				                                        .sub(blockStart)));
			mutableBlocks.get(blockStart)
			             .add(entry.getValue(),
			                  entry.getKey()
			                       .sub(blockStart));
		}
		blocks = ImmutableMap.copyOf(mutableBlocks);
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
		final IntCoordinates blockCoords = getBlock(coordinates);
		final SparseMatrix block = blocks.get(blockCoords.concatenate(blockStarts));
		if (block == null)
			return 0;
		else
			return block.at(coordinates[0] - blockStarts[blockCoords.get(0)],
			                coordinates[1] - blockStarts[blockCoords.get(1)]);
	}
	
	private IntCoordinates getBlock(final int... coordinates)
	{
		
		int blockY = 0;
		int blockX = 0;
		for (int i = 0; i < blockStarts.length; i++)
		{
			if (coordinates[0] >= blockStarts[i])
				blockY = i;
			if (coordinates[1] >= blockStarts[i])
				blockX = i;
		}
		return new IntCoordinates(blockY, blockX);
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
		{
			if (DoubleCompare.almostEqualAfterOps(0,
			                                      this.sub(other)
			                                          .absMaxElement(),
			                                      absmax,
			                                      this.size()))
				return true;
			System.out.println("matrices have different numbers of values. Max difference" + this.sub(other)
			                                                                                     .absMaxElement());
		}
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
		if (PerformanceArguments.getInstance().executeChecks)
			if (getCols() != (vector.getLength()))
				throw new IllegalArgumentException("Incompatible sizes");
		final Map<Integer, Vector> subVectors =
			IntStream.range(0, blockStarts.length)
			         .parallel()
			         .mapToObj(i -> new Tuple2<Integer, Vector>(blockStarts[i],
			                                                    vector.slice(blockStarts[i], blockEnds[i])))
			         .collect(new MapCollector<>());
		final DenseVector ret = new DenseVector(vector.getLength());
		final Map<Integer, List<DenseVector>> multipliedSubVectors =
			blocks.entrySet()
			      .stream()
			      .parallel()
			      .map(e -> new Tuple2<>(e.getKey()
			                              .get(0), e.getValue()
			                                        .mvMul(subVectors.get(
				                                        e.getKey()
				                                         .get(1)))))
			      .collect(Collectors.groupingBy(Tuple2::_1, Collectors.mapping(Tuple2::_2,
			                                                                    Collectors.toList())));
		multipliedSubVectors.forEach((pos, list) -> list.forEach(vec -> ret.addSmallVectorAt(vec, pos)));
		return ret;
	}
	
	@Override
	public Vector tvMul(final Vector vector)
	{
		if (PerformanceArguments.getInstance().executeChecks)
			if (getCols() != (vector.getLength()))
				throw new IllegalArgumentException("Incompatible sizes");
		final Map<Integer, Vector> subVectors =
			IntStream.range(0, blockStarts.length)
			         .parallel()
			         .mapToObj(i -> new Tuple2<Integer, Vector>(blockStarts[i],
			                                                    vector.slice(blockStarts[i], blockEnds[i])))
			         .collect(new MapCollector<>());
		final DenseVector ret = new DenseVector(vector.getLength());
		final Map<Integer, List<DenseVector>> multipliedSubVectors =
			blocks.entrySet()
			      .stream()
			      .parallel()
			      .map(e -> new Tuple2<>(e.getKey()
			                              .get(1), e.getValue()
			                                        .tvMul(subVectors.get(
				                                        e.getKey()
				                                         .get(0)))))
			      .collect(Collectors.groupingBy(Tuple2::_1, Collectors.mapping(Tuple2::_2,
			                                                                    Collectors.toList())));
		multipliedSubVectors.forEach((pos, list) -> list.forEach(vec -> ret.addSmallVectorAt(vec, pos)));
		return ret;
	}
	
	@Override
	public Matrix mmMul(final Matrix matrix)
	{
		return toSparse().mmMul(matrix);
	}
	
	public DenseMatrix mmMul(final DenseMatrix matrix)
	{
		return toSparse().mmMul(matrix);
	}
	
	public SparseMatrix mmMul(final SparseMatrix matrix)
	{
		return toSparse().mmMul(matrix);
	}
	
	public SparseMatrix mmMul(final BlockSparseMatrix matrix)
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
