package dlm;

import basic.CellIntegral;
import basic.RightHandSideIntegral;
import basic.ScalarFunction;
import basic.VectorFunctionOnCells;
import linalg.CoordinateVector;
import linalg.IntCoordinates;
import mixed.MixedCellIntegral;
import mixed.MixedTPCellIntegral;
import mixed.QkQkFunction;
import tensorproduct.ContinuousTPShapeFunction;
import tensorproduct.ContinuousTPVectorFunction;
import tensorproduct.TPVectorCellIntegral;
import tensorproduct.TPVectorCellIntegralOnCell;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

import java.util.ArrayList;
import java.util.List;
import java.util.function.Function;
import java.util.function.Predicate;

public class BackgroundFluidMG
	extends MultiGridFluid
{
	private final double density;
	private final double reynold;
	
	public BackgroundFluidMG(final CoordinateVector startCoordinates,
	                         final CoordinateVector endCoordinates,
	                         final IntCoordinates cells,
	                         final int polynomialDegree,
	                         final int refinements,
	                         final double dt,
	                         final double density,
	                         final double renold)
	{
		super(startCoordinates, endCoordinates, cells, refinements, polynomialDegree);
		this.density = density;
		this.reynold = renold;
	}
	
	@Override
	public List<CellIntegral<TPCell, QkQkFunction>> getIntegrals()
	{
		final TPVectorCellIntegral<ContinuousTPVectorFunction> symGrad =
			new TPVectorCellIntegral<>(1. / reynold, TPVectorCellIntegral.SYM_GRAD);
		final MixedCellIntegral<TPCell, ContinuousTPShapeFunction, ContinuousTPVectorFunction, QkQkFunction>
			symGradMixed = MixedTPCellIntegral.fromVelocityIntegral(symGrad);
		final MixedCellIntegral<TPCell, ContinuousTPShapeFunction, ContinuousTPVectorFunction, QkQkFunction>
			divValue
			= new MixedTPCellIntegral<>(ScalarFunction.constantFunction(-1), MixedTPCellIntegral.DIV_VALUE);
		return List.of(symGradMixed, divValue);
	}
	
	@Override
	public List<CellIntegral<TPCell, QkQkFunction>> getMassIntegrals()
	{
		final TPVectorCellIntegral<ContinuousTPVectorFunction> mass =
			new TPVectorCellIntegral<>(density, TPVectorCellIntegral.VALUE_VALUE);
		final MixedCellIntegral<TPCell, ContinuousTPShapeFunction, ContinuousTPVectorFunction, QkQkFunction>
			mixedMass
			= MixedCellIntegral.fromVelocityIntegral(mass);
		return List.of(mixedMass);
	}
	
	@Override
	public List<CellIntegral<TPCell, QkQkFunction>> getSemiImplicitIntegrals(final VectorFunctionOnCells<TPCell, TPFace> velocity)
	{
		final VectorFunctionOnCells<TPCell, TPFace> semiImplicitWeight1 =
			VectorFunctionOnCells.fromLambda((x) -> velocity.value(x)
			                                                .mul(1. / 2),
			                                 (x, cell) -> velocity.valueInCell(x, cell)
			                                                      .mul(1. / 2), 2, 2);
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
	
	@Override
	public List<RightHandSideIntegral<TPCell, QkQkFunction>> getForceIntegrals(final double t)
	{
		return new ArrayList<>();
	}
	
	@Override
	public Function<CoordinateVector, CoordinateVector> getInitialVelocity()
	{
		return x ->
			CoordinateVector.getUnitVector(2, 0)
			                .mul(0.0);
	}
	
	@Override
	public Function<CoordinateVector, CoordinateVector> velocityBoundaryValues(final double t)
	{
		return x ->
		{
			if (x.x() == 0)
				return CoordinateVector.fromValues(1, 0);
			else
				return CoordinateVector.fromValues(Math.max(0,
				                                            2 - Math.pow(x.dist(CoordinateVector.fromValues(
					                                            0,
					                                            0.5)), 2) * 4), 0);
		};
	}
	
	@Override
	public Predicate<TPFace> getDirichletBoundary()
	{
		return face ->// true;
			face.center()
			    .x() == 0 || face.center()
			                     .y() == 0 || face.center()
			                                      .y() == 1;
	}
}
