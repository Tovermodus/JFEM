package distorted;

import basic.Cell;
import basic.Face;
import basic.FunctionOnCells;
import linalg.CoordinateVector;

public interface ReferenceCellFunction<CT extends Cell<CT, FT>, FT extends Face<CT, FT>, valueT, gradientT, hessianT>
	extends FunctionOnCells<CT, FT, valueT, gradientT, hessianT>
{
	
	@Override
	default valueT valueInCell(final CoordinateVector pos, final CT cell)
	{
		if (cell == null)
			return defaultValue();
		return valueOnReferenceCell(cell.transformToReferenceCell(pos), cell);
	}
	
	@Override
	default gradientT gradientInCell(final CoordinateVector pos, final CT cell)
	{
		if (cell == null)
			return defaultGradient();
		return gradientOnReferenceCell(cell.transformToReferenceCell(pos), cell);
	}
	
	valueT valueOnReferenceCell(final CoordinateVector pos, CT cell);
	
	gradientT gradientOnReferenceCell(final CoordinateVector pos, CT cell);
}
