package basic;

import linalg.CoordinateVector;

public interface FunctionOnCells<CT extends Cell<CT, FT>, FT extends Face<CT, FT>, valueT, gradientT, hessianT>
	extends Function<valueT, gradientT, hessianT>
{
	valueT valueInCell(CoordinateVector pos, CT cell);
	
	gradientT gradientInCell(CoordinateVector pos, CT cell);
	
	default hessianT hessianInCell(final CoordinateVector pos, final CT cell)
	{
		throw new UnsupportedOperationException();
	}
}
