package distorted;

import basic.FunctionOnCells;
import basic.ScalarFunctionOnCells;
import distorted.geometry.DistortedCell;
import distorted.geometry.DistortedFace;
import linalg.CoordinateVector;

public interface DistortedFunction<valueT, gradientT, hessianT>
	extends FunctionOnCells<DistortedCell, DistortedFace, valueT, gradientT, hessianT>
{
	
	@Override
	default valueT valueInCell(final CoordinateVector pos, final DistortedCell cell)
	{
		return valueOnReferenceCell(cell.transformToReferenceCell(pos), cell);
	}
	
	@Override
	default gradientT gradientInCell(final CoordinateVector pos, final DistortedCell cell)
	{
		return gradientOnReferenceCell(cell.transformToReferenceCell(pos), cell);
	}
	
	valueT valueOnReferenceCell(final CoordinateVector pos, DistortedCell cell);
	
	gradientT gradientOnReferenceCell(final CoordinateVector pos, DistortedCell cell);
}
