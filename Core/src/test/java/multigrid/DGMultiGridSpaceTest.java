package multigrid;

import basic.ScalarFESpaceFunction;
import examples.ConvergenceOrderEstimator;
import io.vavr.Tuple2;
import linalg.*;
import org.junit.Test;
import tensorproduct.TPFESpace;
import tensorproduct.TPShapeFunction;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

import java.util.ArrayList;
import java.util.List;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

public class DGMultiGridSpaceTest
{
	
	@Test
	public void testProlongateInterpolate()
	{
		final MGSpace<TPFESpace, TPCell, TPFace, TPShapeFunction, Double,
			CoordinateVector,
			CoordinateMatrix>
			mg = new MGSpace<>(2,
			                   1)
		{
			@Override
			public List<TPFESpace> createSpaces(final int refinements)
			{
				final ArrayList<TPFESpace> ret = new ArrayList<>();
				int mul = 1;
				for (int i = 0; i < refinements + 1; i++)
				{
					ret.add(new TPFESpace(CoordinateVector.fromValues(100.7, .2),
					                      CoordinateVector.fromValues(1000,
					                                                  0.21),
					                      new IntCoordinates(4, 4).mul(mul)));
					mul *= 2;
				}
				return ret;
			}
			
			@Override
			public Tuple2<VectorMultiplyable, DenseVector> createSystem(final TPFESpace space)
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
		final ScalarFESpaceFunction<TPShapeFunction> coars =
			new ScalarFESpaceFunction<>(mg.spaces.get(0)
			                                     .getShapeFunctionMap(), v);
		final ScalarFESpaceFunction<TPShapeFunction> fin =
			new ScalarFESpaceFunction<>(mg.spaces.get(1)
			                                     .getShapeFunctionMap(), prolongator.mvMul(v));
		assertTrue(ConvergenceOrderEstimator.normL2Difference(coars, fin,
		                                                      mg.spaces.get(0)
		                                                               .generatePlotPoints(30)) <= 1e-8);
	}
	
	@Test
	public void testProlongateInterpolateAlgebraic()
	{
		final AMGSpace<TPFESpace, TPCell, TPFace, TPShapeFunction, Double, CoordinateVector,
			CoordinateMatrix>
			mg = new AMGSpace<>(1,
			                    3)
		{
			@Override
			public List<TPFESpace> createSpaces(final int refinements)
			{
				final ArrayList<TPFESpace> ret = new ArrayList<>();
				int mul = 1;
				for (int i = 0; i < refinements + 1; i++)
				{
					ret.add(new TPFESpace(CoordinateVector.fromValues(100.7, .2),
					                      CoordinateVector.fromValues(1000,
					                                                  0.21),
					                      new IntCoordinates(4, 4).mul(mul)));
					mul *= 2;
				}
				return ret;
			}
			
			@Override
			public Tuple2<SparseMatrix, Vector> createFinestLevelSystem(final TPFESpace space)
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
		final ScalarFESpaceFunction<TPShapeFunction> coars =
			new ScalarFESpaceFunction<>(mg.spaces.get(0)
			                                     .getShapeFunctionMap(), v);
		final ScalarFESpaceFunction<TPShapeFunction> fin =
			new ScalarFESpaceFunction<>(mg.spaces.get(1)
			                                     .getShapeFunctionMap(), prolongator.mvMul(v));
		final ScalarFESpaceFunction<TPShapeFunction> coafin =
			new ScalarFESpaceFunction<>(mg.spaces.get(0)
			                                     .getShapeFunctionMap(),
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
