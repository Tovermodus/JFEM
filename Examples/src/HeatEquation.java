import basic.*;
import com.google.common.primitives.Ints;
import linalg.*;
import tensorproduct.*;

import java.util.List;
import java.util.Map;

public class HeatEquation
{
	public static void main(final String[] args)
	{
		final PerformanceArguments.PerformanceArgumentBuilder builder =
			new PerformanceArguments.PerformanceArgumentBuilder();
		builder.build();
		final CoordinateVector start = CoordinateVector.fromValues(-1, -1);
		final CoordinateVector end = CoordinateVector.fromValues(1, 1);
		final int polynomialDegree = 2;
		final TPFESpace grid = new TPFESpace(start, end,
		                                     Ints.asList(25, 25));
		grid.assembleCells();
		grid.assembleFunctions(polynomialDegree);
		final TPCellIntegral<TPShapeFunction> gg =
			new TPCellIntegralViaReferenceCell<TPShapeFunction>(1,
			                                                    TPCellIntegral.GRAD_GRAD);
		final TPCellIntegral<TPShapeFunction> vv =
			new TPCellIntegralViaReferenceCell<TPShapeFunction>(1,
			                                                    TPCellIntegral.VALUE_VALUE);
		final double penalty = 200000;
		final TPFaceIntegral<TPShapeFunction> jj = new TPFaceIntegral<>(
			ScalarFunction.constantFunction(penalty),
			TPFaceIntegral.VALUE_JUMP_VALUE_JUMP);
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
				return 1. / (1 + 1e3 * pos.sub(c).euclidianNorm() * pos.sub(c2).euclidianNorm());
			}
		};
		final TPRightHandSideIntegral<TPShapeFunction> src = new TPRightHandSideIntegral<>(sourceFun
			, TPRightHandSideIntegral.VALUE);
		final double dt = 0.02;
		final int n = grid.getShapeFunctions().size();
		Vector iterate = new DenseVector(n);
		ScalarFESpaceFunction<TPShapeFunction> u_t;
		final SparseMatrix M = new SparseMatrix(n, n);
		grid.writeCellIntegralsToMatrix(List.of(vv), M);
		grid.writeFaceIntegralsToMatrix(List.of(jj), M);
		final SparseMatrix A = new SparseMatrix(n, n);
		grid.writeCellIntegralsToMatrix(List.of(gg), A);
		//grid.writeFaceIntegralsToMatrix(List.of(jj), A);
		A.mulInPlace(dt);
		final Matrix MA = M.add(A);
		final DenseVector source = new DenseVector(n);
		grid.writeCellIntegralsToRhs(List.of(src), source);
		source.mulInPlace(dt);
		System.out.println("source" + source);
		//System.out.println("M"+M);
		//System.out.println("A"+A);
		final int timesteps = 50;
		final Map<CoordinateVector, Double> vals;
		final List<CoordinateVector> points = grid.generatePlotPoints(60);
		vals = (new ScalarFESpaceFunction<>(
			grid.getShapeFunctions(), iterate)
			        .valuesInPointsAtTime(points, 0));
		System.out.println("x" + iterate);
		final boolean implicit = true;
		for (int i = 0; i < timesteps; i++)
		{
			final IterativeSolver its = new IterativeSolver();
			its.showProgress = false;
			if (implicit)
			{
				final DenseVector rhs = source.add(M.mvMul(iterate));
				iterate = its.solveCG(MA, rhs, 1e-10);
			} else
			{
				final DenseVector rhs = source.add(M.mvMul(iterate)).sub(A.mvMul(iterate));
				iterate = its.solveCG(M, rhs, 1e-10);
			}
			//System.out.println("x"+iterate);
			System.out.println(i);
			vals.putAll(new ScalarFESpaceFunction<>(
				grid.getShapeFunctions(), iterate)
				            .valuesInPointsAtTime(points, dt * i));
		}
		final PlotWindow p = new PlotWindow();
		p.addPlot(new ScalarPlot2DTime(vals, 60, ""));
	}
}
