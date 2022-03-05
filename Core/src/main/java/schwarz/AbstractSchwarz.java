package schwarz;

import basic.Cell;
import basic.Face;
import linalg.Vector;
import linalg.VectorMultiplyable;

import java.util.Collection;

public interface AbstractSchwarz<CT extends Cell<CT, FT>, FT extends Face<CT, FT>, OT extends VectorMultiplyable>
{
	
	Collection<CT> getCellPatch(int patch);
	
	OT getGlobalOperator();
	
	OT getRestrictionOperator(int patch);
	
	OT getLocalOperator(int patch);
	
	default Vector getLocalVector(final int patch, final Vector globalVector)
	{
		return getRestrictionOperator(patch).mvMul(globalVector);
	}
	
	default Vector getGlobalVector(final int patch, final Vector globalVector)
	{
		return getRestrictionOperator(patch).tvMul(globalVector);
	}
	
	int getPatchCount();
	
	SubspaceCorrection<OT> getSubspaceCorrection();
	
	Vector solveLocalSystem(int patch, Vector localVector);
}
