import basic.PlotWindow;
import basic.ScalarFunction;
import basic.ScalarPlot2D;
import dlm.BSSmoother2;
import io.vavr.Tuple2;
import linalg.*;
import mixed.*;
import multigrid.AMGPreconditionerSpace;
import multigrid.Smoother;
import tensorproduct.ContinuousTPShapeFunction;
import tensorproduct.ContinuousTPVectorFunction;
import tensorproduct.TPVectorCellIntegral;
import tensorproduct.TPVectorRightHandSideIntegral;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

import java.util.ArrayList;
import java.util.List;

public class StokesAMGBS2
{
	public static void main(final String[] args)
	{
		final int refinements = 3;
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
		final AMGPreconditionerSpace<TaylorHoodSpace, TPCell, TPFace, QkQkFunction,
			MixedValue, MixedGradient, MixedHessian>
			mg = new AMGPreconditionerSpace<>(refinements,
			                                  1)
		{
			@Override
			public List<TaylorHoodSpace> createSpaces(final int refinements)
			{
				final ArrayList<TaylorHoodSpace> ret = new ArrayList<>();
				int mul = 1;
				for (int i = 0; i < refinements + 1; i++)
				{
					final TaylorHoodSpace s = new TaylorHoodSpace(CoordinateVector.fromValues(0, 0),
					                                              CoordinateVector.fromValues(1, 1),
					                                              new IntCoordinates(4,
					                                                                 4).mul(mul));
					ret.add(s);
					mul *= 2;
				}
				return ret;
			}
			
			@Override
			public Tuple2<SparseMatrix, DenseVector> createSystem(final TaylorHoodSpace space)
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
					ret.add(new BSSmoother2(3,
					                        1,
					                        spaces.get(i)
					                              .getVelocitySize()));//, d.getInvertedDiagonalMatrix()));
				}
				return ret;
			}
			
			@Override
			public void applyZeroBoundaryConditions(final TaylorHoodSpace space, final MutableVector vector)
			{
				space.projectOntoBoundaryValues(
					new ComposedMixedFunction(ScalarFunction.constantFunction(0)
					                                        .makeIsotropicVectorFunction()),
					vector);
				vector.set(0, space.getVelocitySize());
			}
		};
		mg.verbose = true;
		
		PlotWindow.addPlot(new ScalarPlot2D(LaplaceReferenceSolution.scalarReferenceSolution(),
		                                    mg.spaces.get(refinements)
		                                             .generatePlotPoints(
			                                             30),
		                                    30, "reference"));
		
		final IterativeSolver it = new IterativeSolver(true);
		it.showProgress = true;
		final Vector solution = it.solvePGMRES(mg.finest_system,
		                                       mg,
		                                       mg.finest_rhs,
		                                       1e-8);
		System.out.println(mg.finest_system.mvMul(solution)
		                                   .sub(mg.finest_rhs)
		                                   .absMaxElement());
		final MixedFunction sol =
			new MixedTPFESpaceFunction<>(
				mg.spaces.get(refinements)
				         .getShapeFunctionMap(), solution);
		PlotWindow.addPlot(new MixedPlot2D(sol,
		                                   mg.spaces.get(refinements)
		                                            .generatePlotPoints(30),
		                                   30));
	}
}
