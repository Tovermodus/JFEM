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

public class BenchmarkFluid
	extends MultiGridFluid
{
	private final double density;
	private final double viscosity;
	private final double um;
	
	public BenchmarkFluid(final CoordinateVector startCoordinates, final CoordinateVector endCoordinates,
	                      final IntCoordinates cells, final int polynomialDegree, final DLMBenchmark.DLMBenchmarkConfig config,
	                      final double density, final double viscosity)
	{
		super(startCoordinates, endCoordinates, cells, config.refinements, polynomialDegree);
		this.um = config.um;
		this.density = density;
		this.viscosity = viscosity;
	}
	
	@Override
	public List<CellIntegral<TPCell, QkQkFunction>> getIntegrals()
	{
		final TPVectorCellIntegral<ContinuousTPVectorFunction> symGrad =
			new TPVectorCellIntegral<>(density * viscosity, TPVectorCellIntegral.GRAD_GRAD);
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
			                                                .mul(density),
			                                 (x, cell) -> velocity.valueInCell(x, cell)
			                                                      .mul(density), 2, 2);
		final TPVectorCellIntegralOnCell<ContinuousTPVectorFunction> convection1 =
			new TPVectorCellIntegralOnCell<>(semiImplicitWeight1, TPVectorCellIntegralOnCell.GRAD_VALUE);
		final MixedCellIntegral<TPCell, ContinuousTPShapeFunction, ContinuousTPVectorFunction,
			QkQkFunction> mixedConvection1 = MixedCellIntegral.fromVelocityIntegral(convection1);
		return List.of(mixedConvection1);
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
		final double charachteristicFunctionT = 1;
		return x ->
		{
			if (x.x() == 0)
				return CoordinateVector.fromValues(charachteristicFunctionT * 4 * um * x.y() * (endCoordinates.y() - x.y()) / (endCoordinates.y() * endCoordinates.y()),
				                                   0);
			else
				return CoordinateVector.fromValues(0, 0);
		};
	}
	
	@Override
	public Predicate<TPFace> getDirichletBoundary()
	{
		return face ->// true;
			face.center()
			    .x() == startCoordinates.x() || face.center()
			                                        .y() == startCoordinates.y() || face.center()
			                                                                            .y() == endCoordinates.y();
	}
}
