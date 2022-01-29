package multigrid;

import basic.ScalarFESpaceFunction;
import basic.ScalarFunction;
import examples.ConvergenceOrderEstimator;
import io.vavr.Tuple2;
import linalg.*;
import org.junit.Test;
import tensorproduct.ContinuousTPFESpace;
import tensorproduct.ContinuousTPShapeFunction;
import tensorproduct.TPCellIntegral;
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
			mg = new MGSpace<>(1,
			                   1)
		{
			@Override
			public List<ContinuousTPFESpace> createSpaces(final int refinements)
			{
				final ArrayList<ContinuousTPFESpace> ret = new ArrayList<>();
				int mul = 1;
				for (int i = 0; i < refinements + 1; i++)
				{
					ret.add(new ContinuousTPFESpace(CoordinateVector.fromValues(100.7, .2),
					                                CoordinateVector.fromValues(1000,
					                                                            0.21),
					                                new IntCoordinates(4, 4).mul(mul)));
					mul *= 2;
				}
				return ret;
			}
			
			@Override
			public Tuple2<VectorMultiplyable, DenseVector> createSystem(final ContinuousTPFESpace space)
			{
				final SparseMatrix s = new SparseMatrix(space.getShapeFunctionMap()
				                                             .size(),
				                                        space.getShapeFunctionMap()
				                                             .size());
				final DenseVector rhs = new DenseVector(space.getShapeFunctionMap()
				                                             .size());
				space.writeCellIntegralsToMatrix(List.of(new TPCellIntegral<>(1,
				                                                              TPCellIntegral.GRAD_GRAD)),
				                                 s);
				space.writeBoundaryValuesTo(ScalarFunction.constantFunction(0),
				                            s,
				                            rhs);
				return new Tuple2<>(s, rhs);
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
			                                     .getShapeFunctionMap(), v);
		final ScalarFESpaceFunction<ContinuousTPShapeFunction> fin =
			new ScalarFESpaceFunction<>(mg.spaces.get(1)
			                                     .getShapeFunctionMap(), prolongator.mvMul(v));
		assertTrue(ConvergenceOrderEstimator.normL2Difference(coars, fin,
		                                                      mg.spaces.get(0)
		                                                               .generatePlotPoints(30)) <= 1e-8);
		
		final DenseVector v1 = new DenseVector(mg.systems.get(0)
		                                                 .getVectorSize());
		final DenseVector v2 = new DenseVector(mg.systems.get(0)
		                                                 .getVectorSize());
		final SparseMatrix P = mg.prolongationMatrices.get(0);
		final SparseMatrix PAP = new SparseMatrix(P.tmMul((SparseMatrix) mg.systems.get(1))
		                                           .mmMul(P));
		int j = 6;
		for (final IntCoordinates c : v.getShape()
		                               .range())
		{
			v1.set(j++, c);
			v2.set(1, c);
		}
		
		mg.spaces.get(0)
		         .forEachBoundaryFace(f ->
		                              {
			                              final var sfs = mg.spaces.get(0)
			                                                       .getFaceSupportMapping()
			                                                       .get(f);
			                              sfs.forEach(fun ->
			                                          {
				                                          if (fun.getNodeFunctional()
				                                                 .usesFace(f))
				                                          {
					                                          v1.set(0,
					                                                 fun.getGlobalIndex());
					                                          v2.set(0,
					                                                 fun.getGlobalIndex());
				                                          }
			                                          });
		                              });
		mg.spaces.get(0)
		         .writeBoundaryValuesTo(ScalarFunction.constantFunction(0),
		                                PAP,
		                                new DenseVector(v1.getLength()));
		final DenseVector Pv1 = P.mvMul(v1);
		final DenseVector Pv2 = P.mvMul(v2);
		assertTrue(PAP.almostEqual((SparseMatrix) mg.systems.get(0)));
		assertTrue(Math.abs(mg.systems.get(1)
		                              .mvMul(Pv1)
		                              .inner(Pv2) - mg.systems.get(0)
		                                                      .mvMul(v1)
		                                                      .inner(v2)) < 1e-8);
		assertTrue(Math.abs(mg.systems.get(1)
		                              .mvMul(P.mvMul(v1))
		                              .inner(P.mvMul(v2)) - mg.systems.get(0)
		                                                              .mvMul(v1)
		                                                              .inner(v2)) < 1e-8);
		assertTrue(Math.abs(PAP.mvMul(v1)
		                       .inner(v2) - mg.systems.get(0)
		                                              .mvMul(v1)
		                                              .inner(v2)) < 1e-8);
	}
	
	@Test
	public void testProlongateInterpolateAlgebraic()
	{
		final AMGSpace<ContinuousTPFESpace, TPCell, TPFace, ContinuousTPShapeFunction, Double, CoordinateVector,
			CoordinateMatrix>
			mg = new AMGSpace<>(1,
			                    1)
		{
			@Override
			public List<ContinuousTPFESpace> createSpaces(final int refinements)
			{
				final ArrayList<ContinuousTPFESpace> ret = new ArrayList<>();
				int mul = 1;
				for (int i = 0; i < refinements + 1; i++)
				{
					ret.add(new ContinuousTPFESpace(CoordinateVector.fromValues(100.7, .2),
					                                CoordinateVector.fromValues(1000,
					                                                            0.21),
					                                new IntCoordinates(4, 4).mul(mul)));
					mul *= 2;
				}
				return ret;
			}
			
			@Override
			public Tuple2<SparseMatrix, Vector> createFinestLevelSystem(final ContinuousTPFESpace space)
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
		final ScalarFESpaceFunction<ContinuousTPShapeFunction> coars =
			new ScalarFESpaceFunction<>(mg.spaces.get(0)
			                                     .getShapeFunctionMap(), v);
		final ScalarFESpaceFunction<ContinuousTPShapeFunction> fin =
			new ScalarFESpaceFunction<>(mg.spaces.get(1)
			                                     .getShapeFunctionMap(), prolongator.mvMul(v));
		final ScalarFESpaceFunction<ContinuousTPShapeFunction> coafin =
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
