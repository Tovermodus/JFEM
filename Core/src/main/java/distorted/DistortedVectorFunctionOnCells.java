package distorted;

import basic.*;
import distorted.geometry.DistortedCell;
import distorted.geometry.DistortedFace;
import linalg.CoordinateMatrix;
import linalg.CoordinateTensor;
import linalg.CoordinateVector;

public interface DistortedVectorFunctionOnCells
	extends DistortedVectorFunction, VectorFunctionOnCells<DistortedCell, DistortedFace>
{
	
	static <CT extends Cell<CT, FT>, FT extends Face<CT, FT>>
	DistortedVectorFunctionOnCells
	concatenate(final VectorFunction outer, final DistortedVectorFunctionOnCells inner)
	{
		if (PerformanceArguments.getInstance().executeChecks)
		{
			if (inner.getRangeDimension() != outer.getDomainDimension())
				throw new IllegalArgumentException("Inner function has the wrong range");
		}
		
		final VectorFunction function = outer.concatenateWith(inner);
		return new DistortedVectorFunctionOnCells()
		{
			@Override
			public CoordinateVector valueOnReferenceCell(final CoordinateVector pos,
			                                             final DistortedCell cell)
			{
				return outer.value(inner.valueOnReferenceCell(pos, cell));
			}
			
			@Override
			public CoordinateMatrix gradientOnReferenceCell(final CoordinateVector pos,
			                                                final DistortedCell cell)
			{
				return inner.gradientOnReferenceCell(pos, cell)
				            .mmMul(outer.gradient(inner.valueOnReferenceCell(pos, cell)));
			}
			
			@Override
			public int getDomainDimension()
			{
				return function.getDomainDimension();
			}
			
			@Override
			public int getRangeDimension()
			{
				return function.getRangeDimension();
			}
			
			@Override
			public CoordinateVector defaultValue()
			{
				return function.defaultValue();
			}
			
			@Override
			public CoordinateMatrix defaultGradient()
			{
				return function.defaultGradient();
			}
			
			@Override
			public CoordinateTensor defaultHessian()
			{
				return function.defaultHessian();
			}
			
			@Override
			public CoordinateVector value(final CoordinateVector pos)
			{
				return function.value(pos);
			}
			
			@Override
			public CoordinateMatrix gradient(final CoordinateVector pos)
			{
				return function.gradient(pos);
			}
			
			@Override
			public VectorFunction concatenateWith(final VectorFunction f)
			{
				return function.concatenateWith(f);
			}
		};
	}
}
