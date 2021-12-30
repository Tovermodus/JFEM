package linalg;

import basic.MapCollector;
import basic.PerformanceArguments;
import com.google.common.collect.ImmutableMap;
import io.vavr.Tuple2;

import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class BlockDenseMatrix
	implements Matrix
{
	public ImmutableMap<IntCoordinates, DenseMatrix> getBlocks()
	{
		return blocks;
	}
	
	public DenseMatrix getBlockMatrix(final IntCoordinates blockNumber)
	{
		return blocks.get(blockNumber.concatenate(blockStarts));
	}
	
	public DenseMatrix getBlockMatrix(final int... blockNumber)
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
	
	public int[] getBlockSizes()
	{
		final int[] sizes = new int[blockStarts.length];
		for (int i = 0; i < sizes.length; i++)
			sizes[i] = blockEnds[i] - blockStarts[i];
		return sizes;
	}
	
	public BlockDenseMatrix subMatrix(final IntCoordinates[] blockNumbers)
	{
		final Map<IntCoordinates, DenseMatrix> newBlocks = new HashMap<>();
		for (int i = 0; i < blockStarts.length; i++)
			newBlocks.put(new IntCoordinates(blockStarts[i], blockStarts[i]),
			              new DenseMatrix(blockEnds[i] - blockStarts[i], blockEnds[i] - blockStarts[i]));
		for (final IntCoordinates blockNumber : blockNumbers)
		{
			newBlocks.put(blockNumber.concatenate(blockStarts), getBlockMatrix(blockNumber));
		}
		return new BlockDenseMatrix(newBlocks, this.shape);
	}
	
	final private ImmutableMap<IntCoordinates, DenseMatrix> blocks;
	final private int[] blockStarts;
	final private int[] blockEnds;
	final IntCoordinates shape;
	
	public BlockDenseMatrix(final Map<IntCoordinates, DenseMatrix> blocks, final int rows, final int cols)
	{
		this(blocks, new IntCoordinates(rows, cols));
	}
	
	public BlockDenseMatrix(final Map<IntCoordinates, DenseMatrix> blocks, final IntCoordinates shape)
	{
		this.shape = shape;
		this.blocks = ImmutableMap.copyOf(blocks);
		final List<Tuple2<Integer, Integer>> blockStartSizes
			= blocks.keySet()
			        .stream()
			        .flatMap((Function<IntCoordinates, Stream<Tuple2<Integer, Integer>>>) coords ->
			        {
				        final List<Tuple2<Integer, Integer>> startSizes = new ArrayList<>();
				        startSizes.add(new Tuple2<>(coords.get(0),
				                                    blocks.get(coords)
				                                          .getRows()));
				        startSizes.add(new Tuple2<>(coords.get(1),
				                                    blocks.get(coords)
				                                          .getCols()));
				        return startSizes.stream();
			        })
			        .distinct()
			        .collect(Collectors.toList());
		checkForBadBlocks(blockStartSizes);
		blockStarts = new int[blockStartSizes.size()];
		blockEnds = new int[blockStartSizes.size()];
		blockStartSizes.sort(Comparator.comparingInt(t -> t._1));
		for (int i = 0; i < blockStartSizes.size(); i++)
		{
			blockStarts[i] = blockStartSizes.get(i)._1;
			blockEnds[i] = blockStartSizes.get(i)._1 + blockStartSizes.get(i)._2;
		}
		blocks.forEach(this::checkIfBlockFits);
	}
	
	private static void checkForBadBlocks(final List<Tuple2<Integer, Integer>> blockStartSizes)
	{
		final Map<Integer, Integer> startSizesMap = new TreeMap<>();
		for (final Tuple2<Integer, Integer> tuple : blockStartSizes)
		{
			if (startSizesMap.containsKey(tuple._1))
				throw new IllegalArgumentException("Blocks do not fit on grid");
			startSizesMap.put(tuple._1, tuple._2);
		}
	}
	
	public BlockDenseMatrix(final BlockSparseMatrix matrix)
	{
		this(matrix.getBlocks()
		           .entrySet()
		           .stream()
		           .map(e -> new Tuple2<>(e.getKey(),
		                                  new DenseMatrix(e.getValue())))
		           .collect(new MapCollector<>()), matrix.shape);
	}
	
	public BlockDenseMatrix(final Matrix s, final int nBlocksPerDir)
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
	
	public BlockDenseMatrix(final Matrix s, final int[] blockStarts)
	{
		this.shape = s.getShape();
		if (PerformanceArguments.getInstance().executeChecks)
		{
			if (s.getRows() != s.getCols())
				throw new IllegalArgumentException("Must be square");
		}
		this.blockStarts = blockStarts;
		blockEnds = new int[blockStarts.length];
		if (blockEnds.length - 1 >= 0) System.arraycopy(blockStarts, 1, blockEnds, 0, blockEnds.length - 1);
		blockEnds[blockEnds.length - 1] = s.getCols();
		final HashMap<IntCoordinates, DenseMatrix> mutableBlocks = new HashMap<>();
		for (int i = 0; i < blockStarts.length; i++)
		{
			mutableBlocks.put(new IntCoordinates(blockStarts[i], blockStarts[i]),
			                  new DenseMatrix(blockEnds[i] - blockStarts[i],
			                                  blockEnds[i] - blockStarts[i]));
		}
		for (final Map.Entry<IntCoordinates, Double> entry : s.getCoordinateEntryList()
		                                                      .entrySet())
		{
			final IntCoordinates block = getBlock(entry.getKey()
			                                           .asArray());
			final IntCoordinates blockStart;
			if (block != null)
				blockStart = block.concatenate(blockStarts);
			else throw new IllegalStateException("first start is not 0");
			if (!mutableBlocks.containsKey(blockStart))
				mutableBlocks.put(blockStart,
				                  new DenseMatrix(block.concatenate(blockEnds)
				                                       .sub(blockStart)));
			mutableBlocks.get(blockStart)
			             .add(entry.getValue(),
			                  entry.getKey()
			                       .sub(blockStart));
		}
		blocks = ImmutableMap.copyOf(mutableBlocks);
	}
	
	public BlockDenseMatrix getInvertedDiagonalMatrix()
	{
		Stream<Map.Entry<IntCoordinates, DenseMatrix>> str = blocks.entrySet()
		                                                           .stream();
		if (PerformanceArguments.getInstance().parallelizeThreads)
			str = str.parallel();
		final BlockDenseMatrix ret = new BlockDenseMatrix(str.filter(e ->
		                                                             {
			                                                             return e.getKey()
			                                                                     .get(0) == e.getKey()
			                                                                                 .get(1);
		                                                             })
		                                                     .map(e ->
		                                                          {
			                                                          try
			                                                          {
				                                                          return new Tuple2<>(e.getKey(),
				                                                                              e.getValue()
				                                                                               .inverse());
			                                                          } catch (final Exception exc)
			                                                          {
				                                                          return new Tuple2<>(e.getKey(),
				                                                                              DenseMatrix.identity(
					                                                                              e.getValue()
					                                                                               .getCols()));
			                                                          }
		                                                          })
		                                                     .collect(new MapCollector<>()), this.shape);
		System.out.println(Arrays.toString(ret.blockStarts));
		return ret;
	}
	
	private void checkIfBlockFits(final IntCoordinates k, final DenseMatrix v)
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
		final DenseMatrix block;
		if (blockCoords != null)
			block = blocks.get(blockCoords.concatenate(blockStarts));
		else block = null;
		if (block == null)
			return 0;
		else
			return block.at(coordinates[0] - blockStarts[blockCoords.get(0)],
			                coordinates[1] - blockStarts[blockCoords.get(1)]);
	}
	
	private IntCoordinates getBlock(final int... coordinates)
	{
		
		int blockY = -1;
		int blockX = -1;
		for (int i = 0; i < blockStarts.length; i++)
		{
			if (coordinates[0] >= blockStarts[i] && blockEnds[i] > coordinates[0])
				blockY = i;
			if (coordinates[1] >= blockStarts[i] && blockEnds[i] > coordinates[1])
				blockX = i;
		}
		if (blockX == -1 || blockY == -1)
			return null;
		return new IntCoordinates(blockY, blockX);
	}
	
	@Override
	public DenseMatrix add(final Tensor other)
	{
		return this.toDense()
		           .add(other);
	}
	
	public DenseMatrix toDense()
	{
		final DenseMatrix ret = new DenseMatrix(getShape().get(0), getShape().get(1));
		blocks.forEach((key, value) ->
			               value.getCoordinateEntryList()
			                    .forEach((coord, val) ->
				                             ret.add(val, coord.add(key))));
		return ret;
	}
	
	@Override
	public BlockDenseMatrix mul(final double scalar)
	{
		final Map<IntCoordinates, DenseMatrix> ret = new HashMap<>();
		for (final Map.Entry<IntCoordinates, DenseMatrix> entry : blocks.entrySet())
			ret.put(entry.getKey(), entry.getValue()
			                             .mul(scalar));
		return new BlockDenseMatrix(ret, getShape());
	}
	
	@Override
	public List<Vector> unfoldDimension(final int dimension)
	{
		return this.toDense()
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
		return false;
	}
	
	@Override
	public IntCoordinates getShape()
	{
		return shape;
	}
	
	@Override
	public boolean almostEqual(final Tensor other)
	{
		if (PerformanceArguments.getInstance().executeChecks)
		{
			if (!getShape().equals(other.getShape()))
				throw new IllegalArgumentException("Tensors are of different size");
		}
		return toDense().almostEqual(this);
	}
	
	@Override
	public BlockDenseMatrix transpose()
	{
		final Map<IntCoordinates, DenseMatrix> ret = new HashMap<>();
		for (final Map.Entry<IntCoordinates, DenseMatrix> entry : blocks.entrySet())
			ret.put(new IntCoordinates(entry.getKey()
			                                .get(1), entry.getKey()
			                                              .get(0)),
			        entry.getValue()
			             .transpose());
		return new BlockDenseMatrix(ret, getShape());
	}
	
	@Override
	public Vector mvMul(final Vector vector)
	{
		if (PerformanceArguments.getInstance().executeChecks)
			if (getCols() != (vector.getLength()))
				throw new IllegalArgumentException("Incompatible sizes");
		final Map<Integer, DenseVector> subVectors =
			BlockSparseMatrix.partitionVectorWithStart(vector, blockStarts, blockEnds);
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
		final Map<Integer, DenseVector> subVectors =
			BlockSparseMatrix.partitionVectorWithStart(vector, blockStarts, blockEnds);
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
	public DenseMatrix mmMul(final Matrix matrix)
	{
		return toDense().mmMul(matrix);
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
		Stream<Map.Entry<IntCoordinates, DenseMatrix>> str = blocks
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
