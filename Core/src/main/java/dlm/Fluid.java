package dlm;

import basic.*;
import io.vavr.Tuple2;
import linalg.CoordinateVector;
import linalg.DenseVector;
import linalg.SparseMatrix;
import mixed.*;
import org.jetbrains.annotations.NotNull;
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
	
	Function<CoordinateVector, CoordinateVector> velocityBoundaryValues();
	
	Predicate<TPFace> getDirichletBoundary();
	
	default int getVelocitySize()
	{
		return (int) getSpace().getShapeFunctionMap()
		                       .values()
		                       .stream()
		                       .filter(ComposeMixedShapeFunction::hasVelocityFunction)
		                       .count();
	}
	
	default int getSystemSize()
	{
		return (int) getSpace().getShapeFunctionMap()
		                       .size();
	}
	
	default int getPressureSize()
	{
		return (int) getSpace().getShapeFunctionMap()
		                       .values()
		                       .stream()
		                       .filter(ComposeMixedShapeFunction::hasPressureFunction)
		                       .count();
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
	
	FluidSystem buildSystem(final double t, final FluidIterate iterate);
	
	default Tuple2<SparseMatrix, DenseVector> getBlockRhs(final FluidSystem fs, final double dt)
	{
		return getBlockRhsForSpace(getSpace(), getVelocitySize(), fs, dt);
	}
	
	@NotNull
	default Tuple2<SparseMatrix, DenseVector> getBlockRhsForSpace(final TaylorHoodSpace space,
	                                                              final int velocitySize,
	                                                              final FluidSystem fs, final double dt)
	{
		final SparseMatrix s =
			new SparseMatrix(fs.massMatrix.mul(1. / dt)
			                              .add(fs.flowMatrix)
			                              .add(fs.semiImplicitMatrix));
		
		final DenseVector d = new DenseVector(fs.forceRhs.add(fs.accelerationRhs.mul(1. / dt)));
		final MixedFunction velocityBoundary =
			new ComposedMixedFunction(VectorFunction.fromLambda(velocityBoundaryValues(), 2, 2));
		space.writeBoundaryValuesTo(velocityBoundary,
		                            getDirichletBoundary(),
		                            //(face, fun) -> fun.hasVelocityFunction(),
		                            s,
		                            d);
		space.overWriteValue(velocitySize, 0, s, d);
		return new Tuple2<>(s, d);
	}
}
