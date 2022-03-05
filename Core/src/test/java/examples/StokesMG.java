package examples;

import basic.ScalarFunction;
import dlm.BSSmoother5;
import io.vavr.Tuple2;
import linalg.*;
import mixed.*;
import multigrid.MGPreconditionerSpace;
import multigrid.Smoother;
import org.junit.Test;
import tensorproduct.ContinuousTPShapeFunction;
import tensorproduct.ContinuousTPVectorFunction;
import tensorproduct.TPVectorCellIntegral;
import tensorproduct.TPVectorRightHandSideIntegral;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

import java.util.ArrayList;
import java.util.List;

import static org.junit.Assert.assertTrue;

public class StokesMG
{
	
	@Test(timeout = 60000)
	public void testConvergenceBoundary()
	{
		
		for (int refinements = 0; refinements <= 3; refinements++)
		{
			final TPVectorCellIntegral<ContinuousTPVectorFunction> gradGrad =
				new TPVectorCellIntegral<>(ScalarFunction.constantFunction(StokesReferenceSolution.reynolds),
				                           TPVectorCellIntegral.GRAD_GRAD);
			final MixedCellIntegral<TPCell, ContinuousTPShapeFunction, ContinuousTPVectorFunction, QkQkFunction>
				divValue =
				new MixedTPCellIntegral<>(ScalarFunction.constantFunction(-1),
				                          MixedTPCellIntegral.DIV_VALUE);
			final MixedCellIntegral<TPCell, ContinuousTPShapeFunction, ContinuousTPVectorFunction,
				QkQkFunction> vv
				= MixedCellIntegral.fromVelocityIntegral(gradGrad);
			final MixedRightHandSideIntegral<TPCell, ContinuousTPShapeFunction,
				ContinuousTPVectorFunction, QkQkFunction> rightHandSideIntegral =
				MixedRightHandSideIntegral.fromVelocityIntegral(
					new TPVectorRightHandSideIntegral<>(StokesReferenceSolution.rightHandSide(),
					                                    TPVectorRightHandSideIntegral.VALUE));
			final MGPreconditionerSpace<TaylorHoodSpace, TPCell, TPFace, QkQkFunction,
				MixedValue, MixedGradient, MixedHessian>
				mg = new MGPreconditionerSpace<>(refinements,
				                                 1)
			{
				
				@Override
				public List<TaylorHoodSpace> createSpaces(final int refinements)
				{
					final ArrayList<TaylorHoodSpace> ret = new ArrayList<>();
					int mul = 1;
					for (int i = 0; i < refinements + 1; i++)
					{
						final TaylorHoodSpace s
							= new TaylorHoodSpace(CoordinateVector.fromValues(0, 0),
							                      CoordinateVector.fromValues(1, 1),
							                      new IntCoordinates(4,
							                                         4).mul(mul));
						ret.add(s);
						mul *= 2;
					}
					return ret;
				}
				
				@Override
				public Tuple2<VectorMultiplyable, DenseVector> createSystem(
					final TaylorHoodSpace space)
				{
					final SparseMatrix s = new SparseMatrix(space.getShapeFunctionMap()
					                                             .size(),
					                                        space.getShapeFunctionMap()
					                                             .size());
					final DenseVector rhs = new DenseVector(space.getShapeFunctionMap()
					                                             .size());
					space.writeCellIntegralsToMatrix(List.of(vv, divValue), s);
					space.writeCellIntegralsToRhs(List.of(rightHandSideIntegral), rhs);
					space.writeBoundaryValuesTo(new ComposedMixedFunction(StokesReferenceSolution.vectorBoundaryValues()),
					                            f -> true,
					                            (f, sf) -> sf.hasVelocityFunction(),
					                            s,
					                            rhs);
					final int firstPressure =
						space.getVelocitySize();
					space.overWriteValue(firstPressure, 0, s, rhs);
					return new Tuple2<>(s, rhs);
				}
				
				@Override
				public List<Smoother> createSmoothers()
				{
					final ArrayList<Smoother> ret = new ArrayList<>();
					for (int i = 1; i < spaces.size(); i++)
					{
						ret.add(new BSSmoother5(3, 6, 1, getSpace(i).getVelocitySize(), this));
					}
					return ret;
				}
				
				@Override
				public void applyZeroBoundaryConditions(final TaylorHoodSpace space,
				                                        final MutableVector vector)
				{
					space.projectOntoBoundaryValues(
						new ComposedMixedFunction(ScalarFunction.constantFunction(0)
						                                        .makeIsotropicVectorFunction()),
						f -> true,
						(f, sf) -> sf.hasVelocityFunction(),
						vector);
					vector.set(0, space.getVelocitySize());
				}
			};
			mg.verbose = false;
			
			Vector iterate = new DenseVector(mg.finest_rhs.getLength());
			for (int i = 0; i < iterate.getLength(); i++)
				((DenseVector) iterate).
					
					set(Math.random(), i);
			mg.applyZeroBoundaryConditions(mg.getFinestSpace(), ((DenseVector) iterate));
			final double initialres = mg.getFinestSystem()
			                            .mvMul(iterate)
			                            .sub(mg.finest_rhs)
			                            .euclidianNorm();
			double res = initialres;
			int i;
			for (i = 0; i < 100 && res > initialres * 1e-8; i++)
			
			{
				iterate = mg.vCycle(iterate, mg.finest_rhs);
				res = mg.finest_system.mvMul(iterate)
				                      .sub(mg.finest_rhs)
				                      .euclidianNorm();
				System.out.println(res / initialres + " after iterations " + i);
			}
			System.out.println(i);
			assertTrue(i < 7);
			final MixedTPFESpaceFunction<QkQkFunction> solution =
				new MixedTPFESpaceFunction<>(mg.getFinestSpace()
				                               .getShapeFunctionMap(), iterate);
			final double error
				= ConvergenceOrderEstimator.normL2VecDifference(solution.getVelocityFunction(),
				                                                StokesReferenceSolution.velocityReferenceSolution(),
				                                                mg.getFinestSpace()
				                                                  .generatePlotPoints(20));
			assertTrue("error " + error + " > " + 0.2 * Math.pow(0.5, refinements),
			           error < 0.2 * Math.pow(0.5, refinements * 2));
		}
	}
}
