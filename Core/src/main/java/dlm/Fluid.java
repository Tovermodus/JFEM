package dlm;

import basic.*;
import io.vavr.Tuple2;
import it.unimi.dsi.fastutil.ints.Int2DoubleArrayMap;
import it.unimi.dsi.fastutil.ints.Int2DoubleMap;
import linalg.CoordinateVector;
import linalg.DenseVector;
import linalg.SparseMatrix;
import mixed.*;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

import java.util.List;
import java.util.function.Function;
import java.util.function.Predicate;

public interface Fluid
{
	TaylorHoodSpace getSpace();
	
	List<CellIntegral<TPCell, QkQkFunction>>
	getIntegrals();
	
	List<CellIntegral<TPCell, QkQkFunction>>
	getMassIntegrals();
	
	List<CellIntegral<TPCell, QkQkFunction>>
	getSemiImplicitIntegrals(final VectorFunctionOnCells<TPCell, TPFace> velocity);
	
	List<RightHandSideIntegral<TPCell, QkQkFunction>> getForceIntegrals(double t);
	
	default MixedFunctionOnCells<TPCell, TPFace> getVelocityPressure(final FluidIterate iterate)
	{
		return new MixedTPFESpaceFunction<>(getSpace().getShapeFunctionMap(),
		                                    iterate.current);
	}
	
	default VectorFunctionOnCells<TPCell, TPFace> getVelocity(final FluidIterate iterate)
	{
		return getVelocityPressure(iterate).getVelocityFunction();
	}
	
	default ScalarFunctionOnCells<TPCell, TPFace> getPressure(final FluidIterate iterate)
	{
		return getVelocityPressure(iterate).getPressureFunction();
	}
	
	Function<CoordinateVector, CoordinateVector> getInitialVelocity();
	
	Function<CoordinateVector, CoordinateVector> velocityBoundaryValues(double t);
	
	Predicate<TPFace> getDirichletBoundary();
	
	default int getSystemSize()
	{
		return (int) getSpace().getShapeFunctions()
		                       .size();
	}
	
	default FluidIterate buildInitialIterate()
	{
		final DenseVector initial = new DenseVector(getSystemSize());
		final MixedFunction initialVelo
			= new ComposedMixedFunction(VectorFunction.fromLambda(getInitialVelocity(),
			                                                      2, 2));
		getSpace().getShapeFunctionMap()
		          .values()
		          .forEach(function ->
		                   {
			                   initial.set(function.getNodeFunctional()
			                                       .evaluate(initialVelo), function.getGlobalIndex());
		                   });
		return new FluidIterate(initial);
	}
	
	FluidSystem buildSystem(final double t, final FluidIterate iterate, List<Particle> particles);
	
	static Tuple2<SparseMatrix, DenseVector> getBlockRhs(final FluidSystem fs, final double dt)
	{
		final SparseMatrix s =
			new SparseMatrix(fs.massMatrix.mul(1. / dt)
			                              .add(fs.flowMatrix)
			                              .add(fs.semiImplicitMatrix));
		
		final DenseVector d = new DenseVector(fs.forceRhs.add(fs.accelerationRhs.mul(1. / dt)));
		return new Tuple2<>(s, d);
	}
	
	default Int2DoubleMap getDirichletNodeValues(final double t)
	{
		return getDirichletNodeValuesForSpace(getSpace(), t);
	}
	
	default Int2DoubleMap getDirichletNodeValuesForSpace(final TaylorHoodSpace space, final double t)
	{
		final var shapeFunctionMap = space.getShapeFunctionMap();
		final ComposedMixedFunction cmf =
			new ComposedMixedFunction(VectorFunction.fromLambda(velocityBoundaryValues(t), 2, 2));
		final Int2DoubleMap ret = new Int2DoubleArrayMap();
		final int[] nodes = space.getBoundaryNodes(getDirichletBoundary(),
		                                           (tpFace, qkQkFunction) -> qkQkFunction.hasVelocityFunction());
		for (final int node : nodes)
		{
			ret.put(node,
			        shapeFunctionMap.get(node)
			                        .getNodeFunctional()
			                        .evaluate(cmf));
		}
		ret.put(space.getVelocitySize(), 0);
		return ret;
	}
}
