import basic.PlotWindow;
import basic.ScalarPlot2D;
import basic.VectorFunction;
import basic.VectorFunctionOnCells;
import com.google.common.primitives.Ints;
import linalg.*;
import tensorproduct.*;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

import java.util.List;
import java.util.function.BiConsumer;

public class Turing
{
	final ContinuousTPFEVectorSpace grid;
	
	IterativeSolver its;
	SparseMatrix M;
	SparseMatrix A;
	
	public Turing(final double dt, final int timeSteps)
	{
//		final CoordinateVector start = CoordinateVector.fromValues(-Math.PI, -Math.PI);
//		final CoordinateVector end = CoordinateVector.fromValues(Math.PI, Math.PI);
		final CoordinateVector start = CoordinateVector.fromValues(0, 0);
		final CoordinateVector end = CoordinateVector.fromValues(3, 3);
		final int polynomialDegree = 1;
		grid = new ContinuousTPFEVectorSpace(start, end,
		                                     Ints.asList(30, 30));
		grid.assembleCells();
		grid.assembleFunctions(polynomialDegree);
		final int n = grid.getShapeFunctionMap()
		                  .size();
		final TPVectorCellIntegral<ContinuousTPVectorFunction> gg =
			new TPVectorCellIntegral<>("TP");
		final TPVectorCellIntegral<ContinuousTPVectorFunction> vv =
			new TPVectorCellIntegral<>(TPVectorCellIntegral.VALUE_VALUE);
		M = new SparseMatrix(n, n);
		A = new SparseMatrix(n, n);
		its = new IterativeSolver();
		grid.writeCellIntegralsToMatrix(List.of(vv), M);
		grid.writeCellIntegralsToMatrix(List.of(gg), A);
		M.mulInPlace(1. / dt);
		final Matrix sys = M.add(A);
		System.out.println("built");
		DenseVector currentIterate = initializeInitialIterate();
		
		PlotWindow.addPlotShow(new ScalarPlot2D(getRFunction(currentIterate).getComponentFunction(0),
		                                        grid.generatePlotPoints(100), 100));
		final IterativeSolver it = new IterativeSolver(true);
		for (int i = 0; i < timeSteps; i++)
		{
			currentIterate =
				new DenseVector(it.solveCG(sys, M.mvMul(currentIterate)
				                                 .add(getRhs(currentIterate)), 1e-8));
			final VectorFunction f = asFunction(currentIterate);
//			PlotWindow.addPlotShow(new ScalarPlot2D(ScalarFunction.fromLambda(x -> Math.pow(f.value(x)
//			                                                                                 .y(),
//			                                                                                2) * f.value(x)
//			                                                                                      .x(), 2),
//			                                        grid.generatePlotPoints(40), 40, " " + i));
			PlotWindow.addPlotShow(new ScalarPlot2D(f.getComponentFunction(1),
			                                        grid.generatePlotPoints(40),
			                                        40,
			                                        " " + i));
		}
	}
	
	public static void main(final String[] args)
	{
		final Turing integrator = new Turing(0.1, 3000);
	}
	
	protected DenseVector initializeInitialIterate()
	{
		final DenseVector ret = new DenseVector(grid.getShapeFunctionMap()
		                                            .size());
		for (final ContinuousTPVectorFunction sf : grid.getShapeFunctions())
		{
			if (sf.getComponent() == 0)
			{
				ret.set(1, sf.getGlobalIndex());
				if (sf.getNodeFunctionalPoint()
				      .sub(CoordinateVector.fromValues(0.5, 0.5))
				      .absMaxElement() < 0.2)
					ret.set(0.6, sf.getGlobalIndex());
			} else
			{
				ret.set(0 + Math.random() * 0.2, sf.getGlobalIndex());
				if (sf.getNodeFunctionalPoint()
				      .sub(CoordinateVector.fromValues(0.5, 0.5))
				      .absMaxElement() < 0.2)
					ret.set(1, sf.getGlobalIndex());
			}
		}
		return ret;
	}
	
	static VectorFunction boundaryValues()
	{
		return VectorFunction.fromLambda(x ->
		                                 {
			                                 if (Math.abs(x.x() - Math.PI) < 1e-8)
				                                 return CoordinateVector.fromValues(Math.cos(x.y() * 2) + 1,
				                                                                    Math.cos(x.y() * 2) + 1)
				                                                        .mul(0.5);
			                                 if (Math.abs(x.x() + Math.PI) < 1e-8)
				                                 return CoordinateVector.fromValues(Math.cos(x.y() * 4) + 1,
				                                                                    Math.cos(x.y() * 4) + 1)
				                                                        .mul(0.5);
			                                 else
				                                 return CoordinateVector.fromValues(Math.cos(x.x() * 4) + 1,
				                                                                    Math.cos(x.x() * 4) + 1)
				                                                        .mul(0.5);
		                                 }
			, 2,
			                         2);
	}
	
	protected static BiConsumer<Matrix, MutableVector> boundaryApplier()
	{
		return (mat, vec) ->
		{
			//grid.writeBoundaryValuesTo(boundaryValues(), (MutableMatrix) mat, vec);
		};
	}
	
	protected DenseVector getRhs(final DenseVector currentIterate)
	{
		final DenseVector ret = new DenseVector(grid.getShapeFunctionMap()
		                                            .size());
		final VectorFunctionOnCells<TPCell, TPFace> Rfunction = getRFunction(currentIterate);
		final TPVectorRightHandSideIntegralOnCell<ContinuousTPVectorFunction> Rrhs =
			new TPVectorRightHandSideIntegralOnCell<>(Rfunction, TPVectorRightHandSideIntegral.VALUE);
		grid.writeCellIntegralsToRhs(List.of(Rrhs), ret);
		return ret;//.sub(A.mvMul(currentIterate));
	}
	
	private VectorFunctionOnCells<TPCell, TPFace> getRFunction(final DenseVector currentIterate)
	{
		final double f = 0.156;
		final double k = 0.142;
		final VectorFunctionOnCells<TPCell, TPFace> valFunction = asFunction(currentIterate);
		return new VectorFunctionOnCells<TPCell, TPFace>()
		{
			@Override
			public CoordinateVector valueInCell(final CoordinateVector pos, final TPCell cell)
			{
				final double u = valFunction.valueInCell(pos, cell)
				                            .x();
				final double v = valFunction.valueInCell(pos, cell)
				                            .y();
//				System.out.println(u + " uv  " + v);
//				System.out.println(valFunction.value(pos) + " uv  " + v);
				return CoordinateVector.fromValues(-u * v * v + f * (1 - u),
				                                   u * v * v - (k) * v * v);
			}
			
			@Override
			public int getRangeDimension()
			{
				return 2;
			}
			
			@Override
			public int getDomainDimension()
			{
				return 2;
			}
			
			@Override
			public CoordinateVector value(final CoordinateVector pos)
			{
				final double u = valFunction.value(pos)
				                            .x();
				final double v = valFunction.value(pos)
				                            .y();
//				System.out.println(u + " uv  " + v);
//				System.out.println(valFunction.value(pos) + " uv  " + v);
				return CoordinateVector.fromValues(-u * v * v + f * (1 - u),
				                                   u * v * v - (f + k) * v);
			}
		};
	}
	
	private VectorTPFESpaceFunction<ContinuousTPVectorFunction> asFunction(final DenseVector d)
	{
		return new VectorTPFESpaceFunction<>(grid.getShapeFunctionMap(), d);
	}
}
