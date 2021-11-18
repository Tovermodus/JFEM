package mixed;

import linalg.Vector;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

import java.util.Map;

public class MixedTPFESpaceFunction<MF extends ComposeMixedShapeFunction<TPCell, TPFace, ?, ?>>
	extends MixedFESpaceFunction<MF, TPCell, TPFace>
{
	public MixedTPFESpaceFunction(final MF[] functions, final double[] coefficients)
	{
		super(functions, coefficients);
	}
	
	public MixedTPFESpaceFunction(final Map<Integer, MF> functions, final Vector coefficients)
	{
		super(functions, coefficients);
	}
}
