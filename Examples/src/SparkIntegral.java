import distorted.DistortedRightHandSideIntegral;
import distorted.DistortedSpace;
import linalg.CoordinateVector;
import linalg.DenseVector;
import linalg.SparseMatrix;
import org.apache.spark.SparkConf;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.apache.spark.api.java.function.Function2;
import org.apache.spark.api.java.function.PairFlatMapFunction;
import org.apache.spark.api.java.function.PairFunction;
import scala.Tuple2;
import tensorproduct.ContinuousTPFESpace;
import tensorproduct.ContinuousTPShapeFunction;

import java.util.ArrayList;
import java.util.List;

public class SparkIntegral
{
	private static Tuple2<ContinuousTPFESpace, DistortedSpace> generateSpaces()
	{
		final ContinuousTPFESpace space =
			new ContinuousTPFESpace(CoordinateVector.fromValues(0, 0),
			                        CoordinateVector.fromValues(1, 1),
			                        List.of(20, 20));
		space.assembleCells();
		space.assembleFunctions(2);
		final DistortedSpace s = new DistortedSpace(
			CoordinateVector.fromValues(0.5, 0.5), 0.5, 2);
		s.assembleCells();
		s.assembleFunctions(2);
		return new Tuple2<>(space, s);
	}
	
	public static void main(final String[] args)
	{
		final SparkConf sparkConf = new SparkConf().setAppName("JavaWordCount");
		sparkConf.setMaster("local[*]");
		final int nPartitions = 11;
		final JavaSparkContext ctx = new JavaSparkContext(sparkConf);
		final List<Integer> intList = new ArrayList<>();
		for (int i = 0; i < nPartitions; i++)
			intList.add(i);
		final JavaRDD<Integer> partitions = ctx.parallelize(intList);
		final JavaPairRDD<Integer, Tuple2<ContinuousTPFESpace, DistortedSpace>> grids =
			partitions.mapToPair(
				new PairFunction<>()
				{
					@Override
					public Tuple2<Integer, Tuple2<ContinuousTPFESpace, DistortedSpace>> call(final Integer integer)
					{
						
						return new Tuple2<>(integer, generateSpaces());
					}
				});
		final JavaPairRDD<Integer, Tuple2<ContinuousTPFESpace, DistortedSpace>> gridsCached = grids.cache();
		final int nCart =
			gridsCached.mapValues(tuple -> tuple._1.getShapeFunctions().size()).first()._2;
		final int nDist =
			gridsCached.mapValues(tuple -> tuple._2.getShapeFunctions().size()).first()._2;
		final SparseMatrix d =
			gridsCached.flatMapToPair(
				(PairFlatMapFunction<Tuple2<Integer, Tuple2<ContinuousTPFESpace, DistortedSpace>>, Integer, DenseVector>) gridPair ->
				{
					final int gridIndex = gridPair._1;
					final Tuple2<ContinuousTPFESpace, DistortedSpace> grids1 = gridPair._2;
					final ContinuousTPFESpace cGrid = grids1._1;
					final DistortedSpace dGrid = grids1._2;
					final int functionsInFirstPartitions =
						((int) (nCart / nPartitions) * nPartitions) / (nPartitions - 1);
					final int start = gridIndex * functionsInFirstPartitions;
					final int end;
					if (gridIndex != nPartitions - 1)
						end = (gridIndex + 1) * functionsInFirstPartitions;
					else
						end = nCart;
					
					System.out.println(
						gridIndex + " " + start + " " + end + " " + nCart + " " + functionsInFirstPartitions);
					final List<Tuple2<Integer, DenseVector>> vectors = new ArrayList<>(
						end - start);
					for (int i = start; i < end; i++)
					{
						final DenseVector vec = new DenseVector(nDist);
						final DistortedRightHandSideIntegral integral
							=
							new DistortedRightHandSideIntegral(
								cGrid.getShapeFunctions().get(i),
								DistortedRightHandSideIntegral.VALUE);
						dGrid.writeCellIntegralsToRhs(List.of(integral), vec);
						vectors.add(
							new Tuple2<>(cGrid.getShapeFunctions().get(i).getGlobalIndex(),
							             vec));
					}
					return vectors.iterator();
				}).aggregate(new SparseMatrix(nDist, nCart),
			                     (Function2<SparseMatrix, Tuple2<Integer, DenseVector>, SparseMatrix>) (sparseMatrix, indexIntegral) ->
			                     {
				                     sparseMatrix.addColumn(indexIntegral._2,
				                                            indexIntegral._1);
				                     return sparseMatrix;
			                     },
			                     (Function2<SparseMatrix, SparseMatrix, SparseMatrix>) (sparseMatrix, sparseMatrix2) ->
			                     {
				                     sparseMatrix.addInPlace(sparseMatrix2);
				                     return sparseMatrix;
			                     });
		final Tuple2<ContinuousTPFESpace, DistortedSpace> spaces = generateSpaces();
		final SparseMatrix mat = new SparseMatrix(nDist, nCart);
		for (final ContinuousTPShapeFunction sf : spaces._1.getShapeFunctions().values())
		{
			final DistortedRightHandSideIntegral integral
				=
				new DistortedRightHandSideIntegral(
					sf,
					DistortedRightHandSideIntegral.VALUE);
			final DenseVector vec = new DenseVector(nDist);
			//if (sf.getGlobalIndex() == 0)
			spaces._2.writeCellIntegralsToRhs(List.of(integral), vec);
			mat.addColumn(vec, sf.getGlobalIndex());
		}
		System.out.println(mat.sub(d).absMaxElement());
	}
}
