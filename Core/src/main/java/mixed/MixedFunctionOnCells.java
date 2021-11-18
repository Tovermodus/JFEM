package mixed;

import basic.Cell;
import basic.Face;
import basic.FunctionOnCells;

public interface MixedFunctionOnCells<CT extends Cell<CT, FT>, FT extends Face<CT, FT>>
	extends MixedFunction, FunctionOnCells<CT, FT, MixedValue, MixedGradient, MixedHessian>
{

}
