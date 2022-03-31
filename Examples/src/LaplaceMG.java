import basic.PlotWindow;
import basic.ScalarFESpaceFunction;
import basic.ScalarFunction;
import basic.ScalarPlot2D;
import io.vavr.Tuple2;
import linalg.*;
import multigrid.MGPreconditionerSpace;
import multigrid.RichardsonSmoother;
import multigrid.Smoother;
import tensorproduct.CTPFESpace;
import tensorproduct.CTPShapeFunction;
import tensorproduct.TPCellIntegral;
import tensorproduct.TPRightHandSideIntegral;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

import java.util.ArrayList;
import java.util.List;

public class LaplaceMG
{
	public static void main(final String[] args)
	{
		
		final int refinements = 1;
		final TPCellIntegral<CTPShapeFunction> gg =
			new TPCellIntegral<>(TPCellIntegral.GRAD_GRAD);
		final TPRightHandSideIntegral<CTPShapeFunction> rightHandSideIntegral =
			new TPRightHandSideIntegral<>(LaplaceReferenceSolution.scalarRightHandSide(),
			                              TPRightHandSideIntegral.VALUE);
		final MGPreconditionerSpace<CTPFESpace, TPCell, TPFace, CTPShapeFunction, Double,
			CoordinateVector,
			CoordinateMatrix>
			mg = new MGPreconditionerSpace<>(refinements,
			                                 2)
		{
			@Override
			public List<CTPFESpace> createSpaces(final int refinements)
			{
				final ArrayList<CTPFESpace> ret = new ArrayList<>();
				int mul = 1;
				for (int i = 0; i < refinements + 1; i++)
				{
					ret.add(new CTPFESpace(CoordinateVector.fromValues(-1, -1),
					                       CoordinateVector.fromValues(1, 1),
					                       new IntCoordinates(2, 2).mul(mul)));
					mul *= 2;
				}
				return ret;
			}
			
			@Override
			public Tuple2<VectorMultiplyable, DenseVector> createSystem(final CTPFESpace space)
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
			public void applyZeroBoundaryConditions(final CTPFESpace space,
			                                        final MutableVector vector)
			{
				space.projectOntoBoundaryValues(ScalarFunction.constantFunction(0), vector);
			}
		};
		mg.verbose = false;
		DenseVector solut = new DenseVector(mg.finest_rhs.mul(0));
		for (final IntCoordinates c : solut.getShape()
		                                   .range())
			solut.add(Math.random() * 0.1, c);
		PlotWindow.addPlot(new ScalarPlot2D(LaplaceReferenceSolution.scalarReferenceSolution(),
		                                    mg.spaces.get(refinements)
		                                             .generatePlotPoints(
			                                             30),
		                                    30, "reference"));
		ScalarFESpaceFunction<CTPShapeFunction> sol =
			new ScalarFESpaceFunction<>(
				mg.spaces.get(refinements)
				         .getShapeFunctionMap(), solut);
		PlotWindow.addPlot(new ScalarPlot2D(sol,
		                                    mg.spaces.get(refinements)
		                                             .generatePlotPoints(30),
		                                    30, "initial"));
		
		final IterativeSolver it = new IterativeSolver();
		it.showProgress = true;
		solut =
			new DenseVector(new GMRES2(1e-8).solve(mg.finest_system, mg,
			                                       mg.finest_rhs));
		final DenseVector solut2 = new
			DenseVector(it.solvePGMRES(mg.finest_system, mg,
			                           mg.finest_rhs,
			                           1e-8));
		System.out.println("sss" + solut2.sub(solut)
		                                 .absMaxElement());
		System.out.println(mg.finest_system.mvMul(solut)
		                                   .sub(mg.finest_rhs)
		                                   .absMaxElement());
		System.out.println(mg.finest_system.mvMul(solut2)
		                                   .sub(mg.finest_rhs)
		                                   .absMaxElement());
		sol =
			new ScalarFESpaceFunction<>(
				mg.spaces.get(refinements)
				         .getShapeFunctionMap(), solut);
		PlotWindow.addPlot(new ScalarPlot2D(sol,
		                                    mg.spaces.get(refinements)
		                                             .generatePlotPoints(30),
		                                    30, "iteration"));
		System.out.println(solut.getLength());
	}
}
