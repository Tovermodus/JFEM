import basic.ScalarFunction;
import io.vavr.Tuple2;
import linalg.*;
import multigrid.MGPreconditionerSpace;
import multigrid.Smoother;
import schwarz.CartesianUpFrontSchwarz;
import schwarz.DirectSolver;
import schwarz.MultiplicativeSubspaceCorrection;
import schwarz.SchwarzSmoother;
import tensorproduct.ContinuousTPFESpace;
import tensorproduct.ContinuousTPShapeFunction;
import tensorproduct.TPCellIntegral;
import tensorproduct.TPRightHandSideIntegral;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

import java.util.ArrayList;
import java.util.List;

public class LaplaceMGSchwarz
{
	public static void main(final String[] args)
	{
		
		final int refinements = 5;
		final TPCellIntegral<ContinuousTPShapeFunction> gg =
			new TPCellIntegral<>(TPCellIntegral.GRAD_GRAD);
		final TPRightHandSideIntegral<ContinuousTPShapeFunction> rightHandSideIntegral =
			new TPRightHandSideIntegral<>(LaplaceReferenceSolution.scalarRightHandSide(),
			                              TPRightHandSideIntegral.VALUE);
		final MGPreconditionerSpace<ContinuousTPFESpace, TPCell, TPFace, ContinuousTPShapeFunction, Double,
			CoordinateVector,
			CoordinateMatrix>
			mg = new MGPreconditionerSpace<>(refinements,
			                                 2)
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
					System.out.println("space done");
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
				System.out.println("integrating");
				space.writeCellIntegralsToMatrix(List.of(gg), s);
				System.out.println("cell done");
				space.writeCellIntegralsToRhs(List.of(rightHandSideIntegral), rhs);
				System.out.println("rhs done");
				space.writeBoundaryValuesTo(LaplaceReferenceSolution.scalarReferenceSolution(),
				                            s,
				                            rhs);
				System.out.println("bdr done");
				return new Tuple2<>(s, rhs);
			}
			
			@Override
			public List<Smoother> createSmoothers()
			{
				final ArrayList<Smoother> ret = new ArrayList<>();
				for (int i = 1; i < spaces.size(); i++)
				{
					System.out.println("smoother " + i);
					final CartesianUpFrontSchwarz<ContinuousTPShapeFunction> schwarz =
						new CartesianUpFrontSchwarz<>((SparseMatrix) getSystem(i),
						                              getSpace(i),
						                              new IntCoordinates(1, 1)
							                              .mul(Math.max(1, (int) Math.pow(2,
							                                                              i - 2))),
						                              2,
						                              new MultiplicativeSubspaceCorrection<>(),
						                              new DirectSolver());
					ret.add(new SchwarzSmoother(1, schwarz));
					System.out.println("smoother " + i + " done");
				}
				return ret;
			}
			
			@Override
			public void applyZeroBoundaryConditions(final ContinuousTPFESpace space,
			                                        final MutableVector vector)
			{
				space.projectOntoBoundaryValues(ScalarFunction.constantFunction(0), vector);
			}
		};
		mg.verbose = true;
		final IterativeSolver it = new IterativeSolver();
		it.showProgress = true;
		final DenseVector solut2 = new
			DenseVector(it.solvePGMRES(mg.finest_system, mg,
			                           mg.finest_rhs,
			                           1e-8));
		System.out.println(mg.finest_system.mvMul(solut2)
		                                   .sub(mg.finest_rhs)
		                                   .absMaxElement());
	}
}
