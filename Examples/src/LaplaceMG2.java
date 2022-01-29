import basic.PlotWindow;
import basic.ScalarFESpaceFunction;
import basic.ScalarFunction;
import basic.ScalarPlot2D;
import io.vavr.Tuple2;
import linalg.*;
import multigrid.MGPreconditionerSpace;
import multigrid.RichardsonSmoother;
import multigrid.Smoother;
import tensorproduct.ContinuousTPFESpace;
import tensorproduct.ContinuousTPShapeFunction;
import tensorproduct.TPCellIntegral;
import tensorproduct.TPRightHandSideIntegral;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

import java.util.ArrayList;
import java.util.List;

public class LaplaceMG2
{
	public static void main(final String[] args)
	{
		
		final int refinements = 1;
		final TPCellIntegral<ContinuousTPShapeFunction> gg =
			new TPCellIntegral<>(TPCellIntegral.GRAD_GRAD);
		final TPRightHandSideIntegral<ContinuousTPShapeFunction> rightHandSideIntegral =
			new TPRightHandSideIntegral<>(LaplaceReferenceSolution.scalarRightHandSide(),
			                              TPRightHandSideIntegral.VALUE);
		final MGPreconditionerSpace<ContinuousTPFESpace, TPCell, TPFace, ContinuousTPShapeFunction,
			Double,
			CoordinateVector,
			CoordinateMatrix>
			mg = new MGPreconditionerSpace<>(refinements,
			                                 1)
		{
			@Override
			public List<ContinuousTPFESpace> createSpaces(final int refinements)
			{
				final ArrayList<ContinuousTPFESpace> ret = new ArrayList<>();
				int mul = 1;
				for (int i = 0; i < refinements + 1; i++)
				{
					ret.add(new ContinuousTPFESpace(CoordinateVector.fromValues(-1, -1),
					                                CoordinateVector.fromValues(1, 1),
					                                new IntCoordinates(2, 2).mul(mul)));
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
				space.writeCellIntegralsToMatrix(List.of(gg), s);
				space.writeCellIntegralsToRhs(List.of(rightHandSideIntegral), rhs);
				space.writeBoundaryValuesTo(LaplaceReferenceSolution.scalarReferenceSolution(),
				                            s,
				                            rhs);
				return new Tuple2<>(s, rhs);
			}
			
			@Override
			public List<Smoother> createSmoothers()
			{
				final ArrayList<Smoother> ret = new ArrayList<>();
				for (int i = 1; i < spaces.size(); i++)
					ret.add(new RichardsonSmoother(0.1, 5));
				return ret;
			}
			
			@Override
			public void applyCorrectBoundaryConditions(final ContinuousTPFESpace space,
			                                           final MutableVector vector)
			{
				space.projectOntoBoundaryValues(LaplaceReferenceSolution.scalarReferenceSolution(),
				                                vector);
			}
			
			@Override
			public void applyZeroBoundaryConditions(final ContinuousTPFESpace space,
			                                        final MutableVector vector)
			{
				space.projectOntoBoundaryValues(ScalarFunction.constantFunction(0), vector);
			}
		};
		
		mg.verbose = false;
		
		PlotWindow.addPlot(new ScalarPlot2D(LaplaceReferenceSolution.scalarReferenceSolution(),
		                                    mg.spaces.get(refinements)
		                                             .generatePlotPoints(
			                                             30),
		                                    30, "reference"));
		
		final GMRES2 gm = new GMRES2(1e-8, true);
		final GMRES2 gm2 = new GMRES2(1e-8, true);
		gm.verbose = false;
		gm2.verbose = false;
		Vector solution = gm.solve(mg.finest_system, mg.finest_rhs);
		
		solution = gm2.solve(mg.finest_system, mg, mg.finest_rhs);
		System.out.println(mg.finest_system.mvMul(solution)
		                                   .sub(mg.finest_rhs)
		                                   .absMaxElement());
		final ScalarFunction sol =
			new ScalarFESpaceFunction<>(
				mg.spaces.get(refinements)
				         .getShapeFunctionMap(), solution);
		PlotWindow.addPlot(new ScalarPlot2D(sol,
		                                    mg.spaces.get(refinements)
		                                             .generatePlotPoints(30),
		                                    30, "iteration"));
	}
}
