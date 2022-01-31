package examples;

import basic.ScalarFESpaceFunction;
import basic.ScalarFunction;
import io.vavr.Tuple2;
import linalg.*;
import multigrid.*;
import org.junit.Test;
import tensorproduct.*;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

import java.util.ArrayList;
import java.util.List;

import static org.junit.Assert.assertTrue;

public class DGLaplaceMG
{
	@Test(timeout = 30000)
	public void testMG()
	{
		for (int refinements = 1; refinements < 4; refinements++)
		{
			final double penalty = 10000;
			final TPCellIntegral<TPShapeFunction> gg =
				new TPCellIntegral<>(TPCellIntegral.GRAD_GRAD);
			final TPFaceIntegral<TPShapeFunction> jj =
				new TPFaceIntegral<>(penalty, TPFaceIntegral.VALUE_JUMP_VALUE_JUMP);
			final TPRightHandSideIntegral<TPShapeFunction> rightHandSideIntegral =
				new TPRightHandSideIntegral<>(LaplaceReferenceSolution.scalarRightHandSide(),
				                              TPRightHandSideIntegral.VALUE);
			final TPBoundaryFaceIntegral<TPShapeFunction> boundaryFaceIntegral =
				new TPBoundaryFaceIntegral<>(LaplaceReferenceSolution.scalarBoundaryValues(penalty),
				                             TPBoundaryFaceIntegral.VALUE);
			final MGSpace<TPFESpace, TPCell, TPFace, TPShapeFunction, Double,
				CoordinateVector, CoordinateMatrix>
				mg = new MGSpace<>(refinements, 1)
			{
				@Override
				public List<TPFESpace> createSpaces(final int refinements)
				{
					final ArrayList<TPFESpace> ret = new ArrayList<>();
					int mul = 1;
					for (int i = 0; i < refinements + 1; i++)
					{
						ret.add(new TPFESpace(CoordinateVector.fromValues(0, 0),
						                      CoordinateVector.fromValues(1, 1),
						                      new IntCoordinates(2, 2).mul(mul)));
						mul *= 2;
					}
					return ret;
				}
				
				@Override
				public Tuple2<VectorMultiplyable, DenseVector> createSystem(final TPFESpace space)
				{
					final SparseMatrix s = new SparseMatrix(space.getShapeFunctionMap()
					                                             .size(),
					                                        space.getShapeFunctionMap()
					                                             .size());
					final DenseVector rhs = new DenseVector(space.getShapeFunctionMap()
					                                             .size());
					space.writeCellIntegralsToMatrix(List.of(gg), s);
					space.writeCellIntegralsToRhs(List.of(rightHandSideIntegral), rhs);
					space.writeFaceIntegralsToMatrix(List.of(jj), s);
					space.writeFaceIntegralsToRhs(List.of(boundaryFaceIntegral), rhs);
					return new Tuple2<>(s, rhs);
				}
				
				@Override
				public List<Smoother> createSmoothers()
				{
					final ArrayList<Smoother> ret = new ArrayList<>();
					for (int i = 1; i < spaces.size(); i++)
						ret.add(new RichardsonSmoother(0.0003, 30));
					return ret;
				}
			};
			DenseVector solut = new DenseVector(mg.finest_rhs.mul(0));
			for (final IntCoordinates c : solut.getShape()
			                                   .range())
				solut.add(Math.random() * 0.1, c);
			final IterativeSolver it = new IterativeSolver(false);
			solut = new DenseVector(it.solvePGMRES(mg.finest_system,
			                                       new MGPReconditioner2(mg),
			                                       mg.finest_rhs,
			                                       1e-8));
			
			final ScalarFESpaceFunction<TPShapeFunction> sol =
				new ScalarFESpaceFunction<>(
					mg.spaces.get(refinements)
					         .getShapeFunctionMap(), solut);
			final double norm = ConvergenceOrderEstimator
				.normL2Difference(sol,
				                  LaplaceReferenceSolution.scalarReferenceSolution(),
				                  mg.spaces.get(0)
				                           .generatePlotPoints(30));
			assertTrue(0.05 * Math.pow(0.25, refinements) > norm);
		}
		try
		{
			Thread.sleep(10000);
		} catch (final InterruptedException e)
		{
			e.printStackTrace();
		}
	}
	
	@Test(timeout = 30000)
	public void testMGPrecond()
	{
		for (int refinements = 1; refinements < 4; refinements++)
		{
			final double penalty = 10000;
			final TPCellIntegral<TPShapeFunction> gg =
				new TPCellIntegral<>(TPCellIntegral.GRAD_GRAD);
			final TPFaceIntegral<TPShapeFunction> jj =
				new TPFaceIntegral<>(penalty, TPFaceIntegral.VALUE_JUMP_VALUE_JUMP);
			final TPRightHandSideIntegral<TPShapeFunction> rightHandSideIntegral =
				new TPRightHandSideIntegral<>(LaplaceReferenceSolution.scalarRightHandSide(),
				                              TPRightHandSideIntegral.VALUE);
			final TPBoundaryFaceIntegral<TPShapeFunction> boundaryFaceIntegral =
				new TPBoundaryFaceIntegral<>(LaplaceReferenceSolution.scalarBoundaryValues(penalty),
				                             TPBoundaryFaceIntegral.VALUE);
			final MGPreconditionerSpace<TPFESpace, TPCell, TPFace, TPShapeFunction, Double,
				CoordinateVector, CoordinateMatrix>
				mg = new MGPreconditionerSpace<>(refinements, 1)
			{
				@Override
				public List<TPFESpace> createSpaces(final int refinements)
				{
					final ArrayList<TPFESpace> ret = new ArrayList<>();
					int mul = 1;
					for (int i = 0; i < refinements + 1; i++)
					{
						ret.add(new TPFESpace(CoordinateVector.fromValues(0, 0),
						                      CoordinateVector.fromValues(1, 1),
						                      new IntCoordinates(2, 2).mul(mul)));
						mul *= 2;
					}
					return ret;
				}
				
				@Override
				public Tuple2<VectorMultiplyable, DenseVector> createSystem(final TPFESpace space)
				{
					final SparseMatrix s = new SparseMatrix(space.getShapeFunctionMap()
					                                             .size(),
					                                        space.getShapeFunctionMap()
					                                             .size());
					final DenseVector rhs = new DenseVector(space.getShapeFunctionMap()
					                                             .size());
					space.writeCellIntegralsToMatrix(List.of(gg), s);
					space.writeCellIntegralsToRhs(List.of(rightHandSideIntegral), rhs);
					space.writeFaceIntegralsToMatrix(List.of(jj), s);
					space.writeFaceIntegralsToRhs(List.of(boundaryFaceIntegral), rhs);
					return new Tuple2<>(s, rhs);
				}
				
				@Override
				public List<Smoother> createSmoothers()
				{
					final ArrayList<Smoother> ret = new ArrayList<>();
					for (int i = 1; i < spaces.size(); i++)
						ret.add(new RichardsonSmoother(0.0003, 30));
					return ret;
				}
				
				@Override
				public void applyZeroBoundaryConditions(final TPFESpace tpfeSpace,
				                                        final MutableVector vector)
				{
					tpfeSpace.projectOntoBoundaryValues(ScalarFunction.constantFunction(0), vector);
				}
			};
			DenseVector solut = new DenseVector(mg.finest_rhs.mul(0));
			for (final IntCoordinates c : solut.getShape()
			                                   .range())
				solut.add(Math.random() * 0.1, c);
			final IterativeSolver it = new IterativeSolver(false);
			solut = new DenseVector(it.solvePGMRES(mg.finest_system, mg,
			                                       mg.finest_rhs,
			                                       1e-8));
			
			final ScalarFESpaceFunction<TPShapeFunction> sol =
				new ScalarFESpaceFunction<>(
					mg.spaces.get(refinements)
					         .getShapeFunctionMap(), solut);
			final double norm = ConvergenceOrderEstimator
				.normL2Difference(sol,
				                  LaplaceReferenceSolution.scalarReferenceSolution(),
				                  mg.spaces.get(0)
				                           .generatePlotPoints(30));
			assertTrue(0.05 * Math.pow(0.25, refinements) > norm);
		}
		try
		{
			Thread.sleep(10000);
		} catch (final InterruptedException e)
		{
			e.printStackTrace();
		}
	}
	
	@Test(timeout = 30000)
	public void testAMG()
	{
		for (int refinements = 1; refinements < 4; refinements++)
		{
			final double penalty = 10000;
			final TPCellIntegral<TPShapeFunction> gg =
				new TPCellIntegral<>(TPCellIntegral.GRAD_GRAD);
			final TPFaceIntegral<TPShapeFunction> jj =
				new TPFaceIntegral<>(penalty, TPFaceIntegral.VALUE_JUMP_VALUE_JUMP);
			final TPRightHandSideIntegral<TPShapeFunction> rightHandSideIntegral =
				new TPRightHandSideIntegral<>(LaplaceReferenceSolution.scalarRightHandSide(),
				                              TPRightHandSideIntegral.VALUE);
			final TPBoundaryFaceIntegral<TPShapeFunction> boundaryFaceIntegral =
				new TPBoundaryFaceIntegral<>(LaplaceReferenceSolution.scalarBoundaryValues(penalty),
				                             TPBoundaryFaceIntegral.VALUE);
			final AMGSpace<TPFESpace, TPCell, TPFace, TPShapeFunction, Double,
				CoordinateVector, CoordinateMatrix>
				mg = new AMGSpace<>(refinements, 1)
			{
				@Override
				public List<TPFESpace> createSpaces(final int refinements)
				{
					final ArrayList<TPFESpace> ret = new ArrayList<>();
					int mul = 1;
					for (int i = 0; i < refinements + 1; i++)
					{
						ret.add(new TPFESpace(CoordinateVector.fromValues(0, 0),
						                      CoordinateVector.fromValues(1, 1),
						                      new IntCoordinates(2, 2).mul(mul)));
						mul *= 2;
					}
					return ret;
				}
				
				@Override
				public Tuple2<SparseMatrix, Vector> createFinestLevelSystem(final TPFESpace space)
				{
					final SparseMatrix s = new SparseMatrix(space.getShapeFunctionMap()
					                                             .size(),
					                                        space.getShapeFunctionMap()
					                                             .size());
					final DenseVector rhs = new DenseVector(space.getShapeFunctionMap()
					                                             .size());
					space.writeCellIntegralsToMatrix(List.of(gg), s);
					space.writeCellIntegralsToRhs(List.of(rightHandSideIntegral), rhs);
					space.writeFaceIntegralsToMatrix(List.of(jj), s);
					space.writeFaceIntegralsToRhs(List.of(boundaryFaceIntegral), rhs);
					return new Tuple2<>(s, rhs);
				}
				
				@Override
				public List<Smoother> createSmoothers()
				{
					final ArrayList<Smoother> ret = new ArrayList<>();
					for (int i = 1; i < spaces.size(); i++)
						ret.add(new RichardsonSmoother(0.0003, 30));
					return ret;
				}
			};
			DenseVector solut = new DenseVector(mg.finest_rhs.mul(0));
			for (final IntCoordinates c : solut.getShape()
			                                   .range())
				solut.add(Math.random() * 0.1, c);
			final IterativeSolver it = new IterativeSolver(false);
			solut = new DenseVector(it.solvePGMRES(mg.finest_system,
			                                       new MGPReconditioner2(mg),
			                                       mg.finest_rhs,
			                                       1e-8));
			
			final ScalarFESpaceFunction<TPShapeFunction> sol =
				new ScalarFESpaceFunction<>(
					mg.spaces.get(refinements)
					         .getShapeFunctionMap(), solut);
			final double norm = ConvergenceOrderEstimator
				.normL2Difference(sol,
				                  LaplaceReferenceSolution.scalarReferenceSolution(),
				                  mg.spaces.get(0)
				                           .generatePlotPoints(30));
			System.out.println(it.iterations);
			assertTrue(0.05 * Math.pow(0.25, refinements) > norm);
		}
	}
}
