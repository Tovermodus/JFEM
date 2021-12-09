import basic.*;
import com.google.common.primitives.Ints;
import linalg.*;
import scala.Function2;
import tensorproduct.*;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.function.BiConsumer;

public class HeatEquationParabolicIntegrator
	extends FullyImplicitParabolicIntegrator
{
	final TPFESpace grid;
	final ScalarFunction sourceFun = new ScalarFunction()
	{
		final CoordinateVector c = CoordinateVector.fromValues(0.5, 0.5);
		final CoordinateVector c2 = CoordinateVector.fromValues(-0.5, -0.5);
		
		@Override
		public int getDomainDimension()
		{
			return 2;
		}
		
		@Override
		public Double value(final CoordinateVector pos)
		{
			return 1. / (1 + 1e3 * pos.sub(c)
			                          .euclidianNorm() * pos.sub(c2)
			                                                .euclidianNorm());
		}
	};
	
	final TPRightHandSideIntegral<TPShapeFunction> src;
	private final List<CoordinateVector> points;
	IterativeSolver its;
	SparseMatrix M;
	SparseMatrix A;
	Map<CoordinateVector, Double> vals;
	
	public HeatEquationParabolicIntegrator(final double dt, final int timeSteps)
	{
		super(dt, timeSteps);
		final CoordinateVector start = CoordinateVector.fromValues(-1, -1);
		final CoordinateVector end = CoordinateVector.fromValues(1, 1);
		final int polynomialDegree = 2;
		grid = new TPFESpace(start, end,
		                     Ints.asList(25, 25));
		points = grid.generatePlotPoints(30);
		grid.assembleCells();
		grid.assembleFunctions(polynomialDegree);
		final int n = grid.getShapeFunctions()
		                  .size();
		final TPCellIntegral<TPShapeFunction> gg =
			new TPCellIntegralViaReferenceCell<>(TPCellIntegral.GRAD_GRAD);
		final TPCellIntegral<TPShapeFunction> vv =
			new TPCellIntegralViaReferenceCell<>(TPCellIntegral.VALUE_VALUE);
		final double penalty = 200000;
		final TPFaceIntegral<TPShapeFunction> jj = new TPFaceIntegral<>(
			ScalarFunction.constantFunction(penalty),
			TPFaceIntegral.VALUE_JUMP_VALUE_JUMP);
		src = new TPRightHandSideIntegral<>(sourceFun
			, TPRightHandSideIntegral.VALUE);
		M = new SparseMatrix(n, n);
		A = new SparseMatrix(n, n);
		its = new IterativeSolver();
		grid.writeCellIntegralsToMatrix(List.of(vv), M);
		grid.writeCellIntegralsToMatrix(List.of(gg), A);
		grid.writeFaceIntegralsToMatrix(List.of(jj), A);
		vals = new HashMap<>();
	}
	
	public static void main(final String[] args)
	{
		final HeatEquationParabolicIntegrator integrator = new HeatEquationParabolicIntegrator(0.02, 50);
		integrator.loop();
		PlotWindow.addPlot(new ScalarPlot2DTime(integrator.vals, 30, "akdsjh"));
	}
	
	@Override
	protected DenseVector initializeInitialIterate()
	{
		final DenseVector ret = new DenseVector(grid.getShapeFunctions()
		                                            .size());
		vals = asFunction(ret).valuesInPointsAtTime(points, time);
		return ret;
	}
	
	@Override
	protected BiConsumer<Matrix, MutableVector> boundaryApplier()
	{
		return (mat, vec) ->
		{
		};
	}
	
	@Override
	protected DenseVector getAdditionalRhs()
	{
		final DenseVector ret = new DenseVector(grid.getShapeFunctions()
		                                            .size());
		grid.writeCellIntegralsToRhs(List.of(src), ret);
		return ret;
	}
	
	@Override
	protected Matrix getSingleDerivativeOperator()
	{
		return M;
	}
	
	@Override
	protected Matrix getNoDerivativeOperator()
	{
		return A;
	}
	
	@Override
	protected Function2<Matrix, DenseVector, DenseVector> getSolver()
	{
		return (A, b) -> (DenseVector) its.solveCG(A, b, 1e-6);
	}
	
	@Override
	protected void postIterationCallback()
	{
		vals.putAll(asFunction(currentIterate).valuesInPointsAtTime(points, time));
	}
	
	private ScalarFESpaceFunction<TPShapeFunction> asFunction(final DenseVector d)
	{
		return new ScalarFESpaceFunction<>(grid.getShapeFunctions(), d);
	}
}
