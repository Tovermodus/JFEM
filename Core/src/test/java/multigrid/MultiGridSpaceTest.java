package multigrid;

import basic.ScalarFESpaceFunction;
import examples.ConvergenceOrderEstimator;
import io.vavr.Tuple2;
import linalg.*;
import org.junit.Test;
import tensorproduct.ContinuousTPFESpace;
import tensorproduct.ContinuousTPShapeFunction;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

import java.util.ArrayList;
import java.util.List;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

public class MultiGridSpaceTest
{
	
	@Test
	public void testProlongateInterpolate()
	{
		final MGSpace<ContinuousTPFESpace, TPCell, TPFace, ContinuousTPShapeFunction, Double,
			CoordinateVector,
			CoordinateMatrix>
			mg = new MGSpace<>(CoordinateVector.fromValues(100.7, .2),
			                   CoordinateVector.fromValues(1000,
			                                               0.21),
			                   new IntCoordinates(4, 4),
			                   1,
			                   1)
		{
			@Override
			public List<ContinuousTPFESpace> createSpaces(final int refinements)
			{
				final ArrayList<ContinuousTPFESpace> ret = new ArrayList<>();
				int mul = 1;
				for (int i = 0; i < refinements + 1; i++)
				{
					ret.add(new ContinuousTPFESpace(startCoordinates, endCoordinates,
					                                coarseCellsPerDimension.mul(mul)));
					mul *= 2;
				}
				return ret;
			}
			
			@Override
			public Tuple2<VectorMultiplyable, DenseVector> createSystem(final ContinuousTPFESpace space)
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
		final ScalarFESpaceFunction<ContinuousTPShapeFunction> coars =
			new ScalarFESpaceFunction<>(mg.spaces.get(0)
			                                     .getShapeFunctions(), v);
		final ScalarFESpaceFunction<ContinuousTPShapeFunction> fin =
			new ScalarFESpaceFunction<>(mg.spaces.get(1)
			                                     .getShapeFunctions(), prolongator.mvMul(v));
		assertTrue(ConvergenceOrderEstimator.normL2Difference(coars, fin,
		                                                      mg.spaces.get(0)
		                                                               .generatePlotPoints(30)) <= 1e-8);
	}
	
	@Test
	public void testProlongateInterpolateAlgebraic()
	{
		final AMGSpace<ContinuousTPFESpace, TPCell, TPFace, ContinuousTPShapeFunction, Double, CoordinateVector,
			CoordinateMatrix>
			mg = new AMGSpace<>(CoordinateVector.fromValues(100.7, .2),
			                    CoordinateVector.fromValues(1000,
			                                                0.21),
			                    new IntCoordinates(4, 4),
			                    1,
			                    1)
		{
			@Override
			public List<ContinuousTPFESpace> createSpaces(final int refinements)
			{
				final ArrayList<ContinuousTPFESpace> ret = new ArrayList<>();
				int mul = 1;
				for (int i = 0; i < refinements + 1; i++)
				{
					ret.add(new ContinuousTPFESpace(startCoordinates, endCoordinates,
					                                coarseCellsPerDimension.mul(mul)));
					mul *= 2;
				}
				return ret;
			}
			
			@Override
			public Tuple2<SparseMatrix, Vector> createFinestLevelSystem(final ContinuousTPFESpace space)
			{
				return new Tuple2<>(new SparseMatrix(space.getShapeFunctions()
				                                          .size(),
				                                     space.getShapeFunctions()
				                                          .size()),
				                    new DenseVector(space.getShapeFunctions()
				                                         .size()));
			}
			
			@Override
			public List<Smoother> createSmoothers()
			{
				return null;
			}
		};
		System.out.println("createc");
		final SparseMvMul prolongator = mg.prolongationOperator.get(0);
		final SparseMvMul restrictor = mg.restrictionOperator.get(0);
		final DenseVector v = new DenseVector(prolongator.getVectorSize());
		for (final IntCoordinates c : v.getShape()
		                               .range())
			v.set(c.sum(), c);
		final ScalarFESpaceFunction<ContinuousTPShapeFunction> coars =
			new ScalarFESpaceFunction<>(mg.spaces.get(0)
			                                     .getShapeFunctions(), v);
		final ScalarFESpaceFunction<ContinuousTPShapeFunction> fin =
			new ScalarFESpaceFunction<>(mg.spaces.get(1)
			                                     .getShapeFunctions(), prolongator.mvMul(v));
		final ScalarFESpaceFunction<ContinuousTPShapeFunction> coafin =
			new ScalarFESpaceFunction<>(mg.spaces.get(0)
			                                     .getShapeFunctions(),
			                            restrictor.mvMul(prolongator.mvMul(v)));
		assertTrue(ConvergenceOrderEstimator.normL2Difference(coars, coafin,
		                                                      mg.spaces.get(0)
		                                                               .generatePlotPoints(30)) <= 1e-8);
		assertTrue(ConvergenceOrderEstimator.normL2Difference(coars, fin,
		                                                      mg.spaces.get(0)
		                                                               .generatePlotPoints(30)) <= 1e-8);
		assertEquals(restrictor.mvMul(prolongator.mvMul(v)), v);
	}
}