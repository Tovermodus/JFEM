import basic.*;
import dlm.BSSmoother2;
import io.vavr.Tuple2;
import linalg.*;
import mixed.*;
import multigrid.MGPreconditionerSpace;
import multigrid.Smoother;
import org.jetbrains.annotations.NotNull;
import tensorproduct.*;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

public class StokesMGTime
{
	private static MGPreconditionerSpace<TaylorHoodSpace, TPCell, TPFace, QkQkFunction, MixedValue, MixedGradient
		, MixedHessian> getMG(final Vector currentIterate, final double dt, final double time)
	{
		final double reynolds = 100;
		final MixedCellIntegral<TPCell, ContinuousTPShapeFunction, ContinuousTPVectorFunction, QkQkFunction>
			divValue =
			new MixedTPCellIntegral<>(ScalarFunction.constantFunction(1),
			                          MixedTPCellIntegral.DIV_VALUE);
		final MixedCellIntegral<TPCell, ContinuousTPShapeFunction, ContinuousTPVectorFunction, QkQkFunction>
			gradGrad =
			MixedCellIntegral.fromVelocityIntegral(new TPVectorCellIntegral<>(
				ScalarFunction.constantFunction(1. / reynolds),
				TPVectorCellIntegral.GRAD_GRAD));
		final MixedCellIntegral<TPCell, ContinuousTPShapeFunction, ContinuousTPVectorFunction, QkQkFunction>
			valueValue =
			MixedCellIntegral.fromVelocityIntegral(new TPVectorCellIntegral<>(
				ScalarFunction.constantFunction(1),
				TPVectorCellIntegral.VALUE_VALUE));
		return new MGPreconditionerSpace<>(4, 1)
		{
			
			@Override
			public List<TaylorHoodSpace> createSpaces(final int refinements)
			{
				verbose = false;
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
			public Tuple2<VectorMultiplyable, DenseVector> createSystem(final TaylorHoodSpace space)
			{
				final int n = space.getShapeFunctionMap()
				                   .size();
				final SparseMatrix mass = new SparseMatrix(n, n);
				space.writeCellIntegralsToMatrix(List.of(valueValue), mass);
				mass.mulInPlace(1. / dt);
				final SparseMatrix flow = new SparseMatrix(n, n);
				final List<CellIntegral<TPCell, QkQkFunction>> flowIntegrals = new ArrayList<>();
				if (currentIterate != null)
					flowIntegrals.addAll(getConvectionIntegrals(generateCurrentFunction(
						restrictToSize(n, currentIterate),
						space).getVelocityFunction()));
				flowIntegrals.add(gradGrad);
				flowIntegrals.add(divValue);
				space.writeCellIntegralsToMatrix(flowIntegrals, flow);
				final DenseVector src = new DenseVector(n);
				space.writeCellIntegralsToRhs(List.of(getSourceIntegral(time)), src);
				
				final SparseMatrix s = mass.add(flow);
				DenseVector fromLast = new DenseVector(n);
				if (currentIterate != null && currentIterate.getLength() == n)
					fromLast = mass.mvMul(currentIterate);
				final DenseVector d = src.add(fromLast);
				
				space.writeBoundaryValuesTo(new ComposedMixedFunction(getBoundaryFunction(time)),
				                            s,
				                            d);
				
				space.overWriteValue(space.getVelocitySize(), 0, s, d);
				return new Tuple2<>(s, d);
			}
			
			@Override
			public List<Smoother> createSmoothers()
			{
				final ArrayList<Smoother> ret = new ArrayList<>();
				for (int i = 1; i < spaces.size(); i++)
				{
					ret.add(new BSSmoother2(7,
					                        0.3,
					                        spaces.get(i)
					                              .getVelocitySize()));//, d.getInvertedDiagonalMatrix()));
				}
				return ret;
			}
			
			@Override
			public void applyZeroBoundaryConditions(final TaylorHoodSpace space, final MutableVector vector)
			{
				space.projectOntoBoundaryValues(new ComposedMixedFunction(ScalarFunction.constantFunction(
					                                                                        0)
				                                                                        .makeIsotropicVectorFunction()),
				                                vector);
				vector.set(0, space.getVelocitySize());
			}
		};
	}
	
	public static VectorFunction getBoundaryFunction(final double t)
	{
		
		return new VectorFunction()
		{
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
				if (Math.abs(pos.x()) <= 1e-10 || Math.abs(pos.x() - 1) <= 1e-10)
					return CoordinateVector.fromValues(t * 15, t * 40);
				///if(Math.abs(pos.x()-1) <= 1e-1)
				//	return CoordinateVector.fromValues(1,1).mul(Math.sin(t)*Math.sin(t)*10*(0.25-
				//	(0.5-pos.y())*(0.5-pos.y())));
				return new CoordinateVector(2);
			}
		};
	}
	
	public static VectorFunction getInitialFunction()
	{
		
		return new VectorFunction()
		{
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
				return new CoordinateVector(2);
			}
		};
	}
	
	private static Vector getInitialIterate(final TaylorHoodSpace space)
	{
		final DenseVector initial = new DenseVector(space.getShapeFunctions()
		                                                 .size());
		final MixedFunction initialVelo
			= new ComposedMixedFunction(getInitialFunction());
		space.getShapeFunctionMap()
		     .values()
		     .forEach(function ->
		              {
			              initial.set(function.getNodeFunctional()
			                                  .evaluate(initialVelo), function.getGlobalIndex());
		              });
		return initial;
	}
	
	public static MixedFunctionOnCells<TPCell, TPFace> generateCurrentFunction(final Vector iterate,
	                                                                           final TaylorHoodSpace grid)
	{
		return new MixedTPFESpaceFunction<>(grid.getShapeFunctionMap(), iterate);
	}
	
	public static void main(final String[] args)
	{
		
		final PerformanceArguments.PerformanceArgumentBuilder builder =
			new PerformanceArguments.PerformanceArgumentBuilder();
		builder.build();
		final double dt = 0.1;
		final int timesteps = 30;
		final int nPoints = 83;
		
		MGPreconditionerSpace<TaylorHoodSpace, TPCell, TPFace, QkQkFunction, MixedValue, MixedGradient
			, MixedHessian> mg = getMG(null, dt, 0);
		
		Vector iterate = getInitialIterate(mg.getFinestSpace());
		
		final List<CoordinateVector> points = mg.getFinestSpace()
		                                        .generatePlotPoints(nPoints);
		
		final Map<CoordinateVector, Double> pvals =
			(generateCurrentFunction(iterate, mg.getFinestSpace()).pressureValuesInPointsAtTime(points, 0));
		final Map<CoordinateVector, Double> divvals =
			(generateCurrentFunction(iterate, mg.getFinestSpace())
				 .getVelocityFunction()
				 .getDivergenceFunction()
				 .valuesInPointsAtTime(points, 0));
		final Map<CoordinateVector, CoordinateVector> vvals =
			(generateCurrentFunction(iterate, mg.getFinestSpace()).velocityValuesInPointsAtTime(points, 0));
		final IterativeSolver its = new IterativeSolver(true);
		
		for (int i = 1; i <= timesteps; i++)
		{
			its.showProgress = true;
			mg = getMG(iterate, dt, i * dt);
			iterate = its.solvePGMRES(mg.getFinestSystem(), mg, mg.finest_rhs, 1e-7);
			System.out.println("Time " + i * dt + " out of " + timesteps * dt + " done");
			pvals.putAll(generateCurrentFunction(iterate, mg.getFinestSpace())
				             .pressureValuesInPointsAtTime(points, dt * i));
			vvals.putAll(generateCurrentFunction(iterate, mg.getFinestSpace())
				             .velocityValuesInPointsAtTime(points, dt * i));
			divvals.putAll(generateCurrentFunction(iterate, mg.getFinestSpace())
				               .getVelocityFunction()
				               .getDivergenceFunction()
				               .valuesInPointsAtTime(points, dt * i));
		}
		PlotWindow.addPlot(new MixedPlot2DTime(pvals, vvals, nPoints));
		PlotWindow.addPlot(new ScalarPlot2DTime(divvals, nPoints, ""));
	}
	
	@NotNull
	private static MixedRightHandSideIntegral<TPCell, ContinuousTPShapeFunction, ContinuousTPVectorFunction,
		QkQkFunction> getSourceIntegral(final double time)
	{
		return MixedRightHandSideIntegral.fromVelocityIntegral(
			new TPVectorRightHandSideIntegral<>(
				new VectorFunction()
				{
					
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
						return CoordinateVector.fromValues(0,
						                                   0 * 1 * Math.sin((pos.x() - time) * 10));
					}
				},
				TPVectorRightHandSideIntegral.VALUE));
	}
	
	public static List<CellIntegral<TPCell, QkQkFunction>> getConvectionIntegrals(final VectorFunctionOnCells<TPCell,
		TPFace> velocity)
	{
		final VectorFunctionOnCells<TPCell, TPFace> semiImplicitWeight1 =
			VectorFunctionOnCells.fromLambda((x) -> velocity.value(x)
			                                                .mul(1 / 2),
			                                 (x, cell) -> velocity.valueInCell(x, cell)
			                                                      .mul(1 / 2), 2, 2);
		final VectorFunctionOnCells<TPCell, TPFace> semiImplicitWeight2 =
			VectorFunctionOnCells.fromLambda((x) -> velocity.value(x)
			                                                .mul(-1. / 2),
			                                 (x, cell) -> velocity.valueInCell(x, cell)
			                                                      .mul(-1. / 2), 2, 2);
		final TPVectorCellIntegralOnCell<ContinuousTPVectorFunction> convection1 =
			new TPVectorCellIntegralOnCell<>(semiImplicitWeight1, TPVectorCellIntegral.GRAD_VALUE);
		final TPVectorCellIntegralOnCell<ContinuousTPVectorFunction> convection2 =
			new TPVectorCellIntegralOnCell<>(semiImplicitWeight2, TPVectorCellIntegral.VALUE_GRAD);
		final MixedCellIntegral<TPCell, ContinuousTPShapeFunction, ContinuousTPVectorFunction,
			QkQkFunction> mixedConvection1 = MixedCellIntegral.fromVelocityIntegral(convection1);
		final MixedCellIntegral<TPCell, ContinuousTPShapeFunction, ContinuousTPVectorFunction,
			QkQkFunction> mixedConvection2 = MixedCellIntegral.fromVelocityIntegral(convection2);
		return List.of(mixedConvection1, mixedConvection2);
	}
}
