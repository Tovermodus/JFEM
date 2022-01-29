package multigrid;

import basic.VectorFESpaceFunction;
import examples.ConvergenceOrderEstimator;
import io.vavr.Tuple2;
import linalg.*;
import org.junit.Test;
import tensorproduct.ContinuousTPFEVectorSpace;
import tensorproduct.ContinuousTPVectorFunction;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

import java.util.ArrayList;
import java.util.List;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

public class VectorMultiGridSpaceTest
{
	@Test
	public void testProlongateInterpolate()
	{
		final MGSpace<ContinuousTPFEVectorSpace, TPCell, TPFace, ContinuousTPVectorFunction, CoordinateVector,
			CoordinateMatrix, CoordinateTensor>
			mg = new MGSpace<>(1, 1)
		{
			@Override
			public List<ContinuousTPFEVectorSpace> createSpaces(final int refinements)
			{
				final ArrayList<ContinuousTPFEVectorSpace> ret = new ArrayList<>();
				int mul = 1;
				for (int i = 0; i < refinements + 1; i++)
				{
					ret.add(new ContinuousTPFEVectorSpace(CoordinateVector.fromValues(100.7, .2),
					                                      CoordinateVector.fromValues(1000,
					                                                                  0.21),
					                                      new IntCoordinates(4, 4).mul(mul)));
					mul *= 2;
				}
				return ret;
			}
			
			@Override
			public Tuple2<VectorMultiplyable, DenseVector> createSystem(final ContinuousTPFEVectorSpace space)
			{
				return new Tuple2<>(null, null);
			}
			
			@Override
			public List<Smoother> createSmoothers()
			{
				return null;
			}
		};
		final SparseMvMul prolongator = mg.prolongationOperators.get(0);
		final DenseVector v = new DenseVector(prolongator.getVectorSize());
		for (final IntCoordinates c : v.getShape()
		                               .range())
			v.set(c.sum(), c);
		final VectorFESpaceFunction<ContinuousTPVectorFunction> coars =
			new VectorFESpaceFunction<>(mg.spaces.get(0)
			                                     .getShapeFunctionMap(), v);
		final VectorFESpaceFunction<ContinuousTPVectorFunction> fin =
			new VectorFESpaceFunction<>(mg.spaces.get(1)
			                                     .getShapeFunctionMap(), prolongator.mvMul(v));
		assertTrue(ConvergenceOrderEstimator.normL2VecDifference(coars, fin,
		                                                         mg.spaces.get(0)
		                                                                  .generatePlotPoints(30)) <= 1e-8);
	}
	
	@Test
	public void testProlongateInterpolateAlgebraic()
	{
		final AMGSpace<ContinuousTPFEVectorSpace, TPCell, TPFace, ContinuousTPVectorFunction, CoordinateVector,
			CoordinateMatrix, CoordinateTensor>
			mg = new AMGSpace<>(1,
			                    2)
		{
			@Override
			public List<ContinuousTPFEVectorSpace> createSpaces(final int refinements)
			{
				final ArrayList<ContinuousTPFEVectorSpace> ret = new ArrayList<>();
				int mul = 1;
				for (int i = 0; i < refinements + 1; i++)
				{
					ret.add(new ContinuousTPFEVectorSpace(CoordinateVector.fromValues(100.7, .2),
					                                      CoordinateVector.fromValues(1000,
					                                                                  0.21),
					                                      new IntCoordinates(4, 4).mul(mul)));
					mul *= 2;
				}
				return ret;
			}
			
			@Override
			public Tuple2<SparseMatrix, Vector> createFinestLevelSystem(final ContinuousTPFEVectorSpace space)
			{
				return new Tuple2<>(new SparseMatrix(space.getShapeFunctionMap()
				                                          .size(),
				                                     space.getShapeFunctionMap()
				                                          .size()),
				                    new DenseVector(space.getShapeFunctionMap()
				                                         .size()));
			}
			
			@Override
			public List<Smoother> createSmoothers()
			{
				return null;
			}
		};
		System.out.println("createc");
		final SparseMvMul prolongator = mg.prolongationOperators.get(0);
		final SparseMvMul restrictor = mg.restrictionOperator.get(0);
		final DenseVector v = new DenseVector(prolongator.getVectorSize());
		for (final IntCoordinates c : v.getShape()
		                               .range())
			v.set(c.sum(), c);
		final VectorFESpaceFunction<ContinuousTPVectorFunction> coars =
			new VectorFESpaceFunction<>(mg.spaces.get(0)
			                                     .getShapeFunctionMap(), v);
		final VectorFESpaceFunction<ContinuousTPVectorFunction> fin =
			new VectorFESpaceFunction<>(mg.spaces.get(1)
			                                     .getShapeFunctionMap(), prolongator.mvMul(v));
		final VectorFESpaceFunction<ContinuousTPVectorFunction> coafin =
			new VectorFESpaceFunction<>(mg.spaces.get(0)
			                                     .getShapeFunctionMap(),
			                            restrictor.mvMul(prolongator.mvMul(v)));
		assertTrue(ConvergenceOrderEstimator.normL2VecDifference(coars, coafin,
		                                                         mg.spaces.get(0)
		                                                                  .generatePlotPoints(30)) <= 1e-8);
		assertTrue(ConvergenceOrderEstimator.normL2VecDifference(coars, fin,
		                                                         mg.spaces.get(0)
		                                                                  .generatePlotPoints(30)) <= 1e-8);
		assertEquals(restrictor.mvMul(prolongator.mvMul(v)), v);
	}
}
