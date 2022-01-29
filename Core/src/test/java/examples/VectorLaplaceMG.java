package examples;

import basic.VectorFESpaceFunction;
import io.vavr.Tuple2;
import linalg.*;
import multigrid.*;
import org.junit.Test;
import tensorproduct.ContinuousTPFEVectorSpace;
import tensorproduct.ContinuousTPVectorFunction;
import tensorproduct.TPVectorCellIntegral;
import tensorproduct.TPVectorRightHandSideIntegral;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

import java.util.ArrayList;
import java.util.List;

import static org.junit.Assert.assertTrue;

public class VectorLaplaceMG
{
	@Test(timeout = 30000)
	public void testMG()
	{
		
		for (int refinements = 2; refinements < 6; refinements++)
		{
			final TPVectorCellIntegral<ContinuousTPVectorFunction> gg =
				new TPVectorCellIntegral<>(TPVectorCellIntegral.GRAD_GRAD);
			final TPVectorRightHandSideIntegral<ContinuousTPVectorFunction> rightHandSideIntegral =
				new TPVectorRightHandSideIntegral<>(LaplaceReferenceSolution.vectorRightHandSide(),
				                                    TPVectorRightHandSideIntegral.VALUE);
			final MGSpace<ContinuousTPFEVectorSpace, TPCell, TPFace, ContinuousTPVectorFunction,
				CoordinateVector, CoordinateMatrix, CoordinateTensor>
				mg = new MGSpace<>(refinements, 1)
			{
				@Override
				public List<ContinuousTPFEVectorSpace> createSpaces(final int refinements)
				{
					final ArrayList<ContinuousTPFEVectorSpace> ret = new ArrayList<>();
					int mul = 1;
					for (int i = 0; i < refinements + 1; i++)
					{
						ret.add(new ContinuousTPFEVectorSpace(CoordinateVector.fromValues(0, 0),
						                                      CoordinateVector.fromValues(1, 1),
						                                      new IntCoordinates(4,
						                                                         4).mul(mul)));
						mul *= 2;
					}
					return ret;
				}
				
				@Override
				public Tuple2<VectorMultiplyable, DenseVector> createSystem(final ContinuousTPFEVectorSpace space)
				{
					final SparseMatrix s = new SparseMatrix(space.getShapeFunctionMap()
					                                             .size(),
					                                        space.getShapeFunctionMap()
					                                             .size());
					final DenseVector rhs = new DenseVector(space.getShapeFunctionMap()
					                                             .size());
					space.writeCellIntegralsToMatrix(List.of(gg), s);
					space.writeCellIntegralsToRhs(List.of(rightHandSideIntegral), rhs);
					space.writeBoundaryValuesTo(LaplaceReferenceSolution.vectorReferenceSolution(),
					                            s,
					                            rhs);
					return new Tuple2<>(s, rhs);
				}
				
				@Override
				public List<Smoother> createSmoothers()
				{
					final ArrayList<Smoother> ret = new ArrayList<>();
					for (int i = 1; i < spaces.size(); i++)
						ret.add(new RichardsonSmoother(0.1, 8));
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
			
			final VectorFESpaceFunction<ContinuousTPVectorFunction> sol =
				new VectorFESpaceFunction<>(
					mg.spaces.get(refinements)
					         .getShapeFunctionMap(), solut);
			assertTrue(it.iterations < 10);
			final double norm = ConvergenceOrderEstimator
				.normL2VecDifference(sol,
				                     LaplaceReferenceSolution.vectorReferenceSolution(),
				                     mg.spaces.get(0)
				                              .generatePlotPoints(30));
			assertTrue(0.2 * Math.pow(0.25, refinements) > norm);
			System.out.println(0.2 * Math.pow(0.25, refinements) + " " + norm + " " + it.iterations);
		}
	}
	
	@Test(timeout = 30000)
	public void testAMG()
	{
		for (int refinements = 2; refinements < 5; refinements++)
		{
			final TPVectorCellIntegral<ContinuousTPVectorFunction> gg =
				new TPVectorCellIntegral<>(TPVectorCellIntegral.GRAD_GRAD);
			final TPVectorRightHandSideIntegral<ContinuousTPVectorFunction> rightHandSideIntegral =
				new TPVectorRightHandSideIntegral<>(LaplaceReferenceSolution.vectorRightHandSide(),
				                                    TPVectorRightHandSideIntegral.VALUE);
			final AMGSpace<ContinuousTPFEVectorSpace, TPCell, TPFace, ContinuousTPVectorFunction,
				CoordinateVector, CoordinateMatrix, CoordinateTensor>
				mg = new AMGSpace<>(refinements, 1)
			{
				@Override
				public List<ContinuousTPFEVectorSpace> createSpaces(final int refinements)
				{
					final ArrayList<ContinuousTPFEVectorSpace> ret = new ArrayList<>();
					int mul = 1;
					for (int i = 0; i < refinements + 1; i++)
					{
						ret.add(new ContinuousTPFEVectorSpace(CoordinateVector.fromValues(0, 0),
						                                      CoordinateVector.fromValues(1, 1),
						                                      new IntCoordinates(4,
						                                                         4).mul(mul)));
						mul *= 2;
					}
					return ret;
				}
				
				@Override
				public Tuple2<SparseMatrix, Vector> createFinestLevelSystem(final ContinuousTPFEVectorSpace space)
				{
					final SparseMatrix s = new SparseMatrix(space.getShapeFunctionMap()
					                                             .size(),
					                                        space.getShapeFunctionMap()
					                                             .size());
					final DenseVector rhs = new DenseVector(space.getShapeFunctionMap()
					                                             .size());
					space.writeCellIntegralsToMatrix(List.of(gg), s);
					space.writeCellIntegralsToRhs(List.of(rightHandSideIntegral), rhs);
					space.writeBoundaryValuesTo(LaplaceReferenceSolution.vectorReferenceSolution(),
					                            s,
					                            rhs);
					return new Tuple2<>(s, rhs);
				}
				
				@Override
				public List<Smoother> createSmoothers()
				{
					final ArrayList<Smoother> ret = new ArrayList<>();
					for (int i = 1; i < spaces.size(); i++)
						ret.add(new RichardsonSmoother(0.1, 4));
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
			
			final VectorFESpaceFunction<ContinuousTPVectorFunction> sol =
				new VectorFESpaceFunction<>(
					mg.spaces.get(refinements)
					         .getShapeFunctionMap(), solut);
			final double norm = ConvergenceOrderEstimator
				.normL2VecDifference(sol,
				                     LaplaceReferenceSolution.vectorReferenceSolution(),
				                     mg.spaces.get(0)
				                              .generatePlotPoints(30));
			assertTrue(0.2 * Math.pow(0.25, refinements) > norm);
			System.out.println(0.2 * Math.pow(0.25, refinements) + " " + norm + " " + it.iterations);
		}
	}
}
