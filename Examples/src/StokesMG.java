import basic.PlotWindow;
import basic.ScalarFunction;
import io.vavr.Tuple2;
import linalg.*;
import mixed.*;
import multigrid.MGSpace;
import multigrid.RichardsonSmoother;
import multigrid.Smoother;
import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import tensorproduct.ContinuousTPShapeFunction;
import tensorproduct.ContinuousTPVectorFunction;
import tensorproduct.TPVectorCellIntegral;
import tensorproduct.TPVectorRightHandSideIntegral;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Optional;

public class StokesMG
{
	public static void main(final String[] args)
	{
		final int refinements = 0;
		final TPVectorCellIntegral<ContinuousTPVectorFunction> gradGrad =
			new TPVectorCellIntegral<>(ScalarFunction.constantFunction(StokesReferenceSolution.reynolds),
			                           TPVectorCellIntegral.GRAD_GRAD);
		final MixedCellIntegral<TPCell, ContinuousTPShapeFunction, ContinuousTPVectorFunction, QkQkFunction>
			divValue =
			new MixedTPCellIntegral<>(ScalarFunction.constantFunction(-1),
			                          MixedTPCellIntegral.DIV_VALUE);
		final MixedCellIntegral<TPCell, ContinuousTPShapeFunction, ContinuousTPVectorFunction, QkQkFunction> vv
			=
			MixedCellIntegral.fromVelocityIntegral(gradGrad);
		final MixedRightHandSideIntegral<TPCell, ContinuousTPShapeFunction,
			ContinuousTPVectorFunction, QkQkFunction> rightHandSideIntegral =
			MixedRightHandSideIntegral.fromVelocityIntegral(
				new TPVectorRightHandSideIntegral<>(StokesReferenceSolution.rightHandSide(),
				                                    TPVectorRightHandSideIntegral.VALUE));
		final MGSpace<TaylorHoodSpace, TPCell, TPFace, QkQkFunction,
			MixedValue, MixedGradient, MixedHessian>
			mg = new MGSpace<>(refinements,
			                   1)
		{
			@Override
			public List<TaylorHoodSpace> createSpaces(final int refinements)
			{
				final ArrayList<TaylorHoodSpace> ret = new ArrayList<>();
				int mul = 1;
				for (int i = 0; i < refinements + 1; i++)
				{
					ret.add(new TaylorHoodSpace(CoordinateVector.fromValues(0, 0),
					                            CoordinateVector.fromValues(1, 1),
					                            new IntCoordinates(4, 4).mul(mul)));
					mul *= 2;
				}
				return ret;
			}
			
			@Override
			public Tuple2<VectorMultiplyable, DenseVector> createSystem(final TaylorHoodSpace space)
			{
				final SparseMatrix s = new SparseMatrix(space.getShapeFunctions()
				                                             .size(),
				                                        space.getShapeFunctions()
				                                             .size());
				final DenseVector rhs = new DenseVector(space.getShapeFunctions()
				                                             .size());
				space.writeCellIntegralsToMatrix(List.of(vv, divValue), s);
				space.writeCellIntegralsToRhs(List.of(rightHandSideIntegral), rhs);
				space.writeBoundaryValuesTo(new ComposedMixedFunction(StokesReferenceSolution.vectorBoundaryValues()),
				                            s,
				                            rhs);
				final Optional<QkQkFunction> firstPressure =
					space.getShapeFunctions()
					     .values()
					     .stream()
					     .filter(ComposeMixedShapeFunction::hasPressureFunction)
					     .findFirst();
				firstPressure.ifPresent(st -> space.overWriteValue(st.getGlobalIndex(),
				                                                   0,
				                                                   s,
				                                                   rhs));
				System.out.println(s.transpose()
				                    .sub(s)
				                    .absMaxElement());
				return new Tuple2<>(s, rhs);
			}
			
			@Override
			public List<Smoother> createSmoothers()
			{
				final ArrayList<Smoother> ret = new ArrayList<>();
				for (int i = 1; i < spaces.size(); i++)
				{
					System.out.println("i");
					final BlockDenseMatrix d = new BlockDenseMatrix((SparseMatrix) systems.get(i),
					                                                spaces.get(i)
					                                                      .getCells()
					                                                      .size() / 20);
					ret.add(new RichardsonSmoother(0.1, 5));//, d.getInvertedDiagonalMatrix()));
				}
				return ret;
			}
		};
		final DenseMatrix A = new DenseMatrix((SparseMatrix) mg.finest_system);
		PlotWindow.addPlot(new MatrixPlot(A));
		PlotWindow.addPlot(new MatrixPlot(A.inverse()));
		System.out.println(A.mul(-1)
		                    .powerIterationSymmetric() * -1);
		final RealMatrix Areal = A.toApacheMatrix();
		final EigenDecomposition ed = new EigenDecomposition(Areal);
		System.out.println(Arrays.toString(ed.getRealEigenvalues()));
		for (int i = 0; i < A.getCols(); i++)
			if (ed.getRealEigenvalue(i) < 0)
			{
				System.out.println(ed.getEigenvector(i));
				final Vector ev = DenseVector.fromRealVector(ed.getEigenvector(i));
				final MixedFESpaceFunction<QkQkFunction, TPCell, TPFace> sol =
					new MixedTPFESpaceFunction<>(
						mg.spaces.get(refinements)
						         .getShapeFunctions(), ev);
				PlotWindow.addPlot(new MixedPlot2D(sol,
				                                   mg.spaces.get(refinements)
				                                            .generatePlotPoints(30),
				                                   30));
			}
//		DenseVector solut = new DenseVector(mg.finest_rhs.mul(0));
//		for (final IntCoordinates c : solut.getShape()
//		                                   .range())
//			solut.add(Math.random() * 0.1, c);
//		PlotWindow.addPlot(new MixedPlot2D(new ComposedMixedFunction(StokesReferenceSolution.pressureReferenceSolution(),
//		                                                             StokesReferenceSolution.velocityReferenceSolution()),
//		                                   mg.spaces.get(refinements)
//		                                            .generatePlotPoints(
//			                                            30),
//		                                   30));
//		MixedFESpaceFunction<QkQkFunction, TPCell, TPFace> sol =
//			new MixedTPFESpaceFunction<>(
//				mg.spaces.get(refinements)
//				         .getShapeFunctions(), solut);
//		PlotWindow.addPlot(new MixedPlot2D(sol,
//		                                   mg.spaces.get(refinements)
//		                                            .generatePlotPoints(30),
//		                                   30));
//
//		final IterativeSolver it = new IterativeSolver();
//		it.showProgress = true;
//
//		solut = new DenseVector(it.solvePCG(new SparseMvMul((SparseMatrix) mg.finest_system),
//		                                    new MGPReconditioner(mg),
//		                                    mg.finest_rhs,
//		                                    1e-5));
//		sol =
//			new MixedTPFESpaceFunction<>(
//				mg.spaces.get(refinements)
//				         .getShapeFunctions(), solut);
//		PlotWindow.addPlot(new MixedPlot2D(sol,
//		                                   mg.spaces.get(refinements)
//		                                            .generatePlotPoints(30),
//		                                   30));
//		System.out.println(solut.getLength());
	}
}
