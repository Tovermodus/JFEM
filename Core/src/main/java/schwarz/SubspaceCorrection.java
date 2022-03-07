package schwarz;

import linalg.Vector;
import linalg.VectorMultiplyable;
import org.jetbrains.annotations.NotNull;

public interface SubspaceCorrection<OT extends VectorMultiplyable>
{
	Vector solve(final AbstractSchwarz<?, ?, OT> schwarz, final Vector globalRhs);
	
	default Vector apply(final AbstractSchwarz<?, ?, OT> schwarz, final @NotNull Vector globalIterate,
	                     final @NotNull Vector globalRhs)
	{
		final Vector globalResidual = globalRhs.sub(schwarz.getGlobalOperator()
		                                                   .mvMul(globalIterate));
		final Vector sol = solve(schwarz, globalResidual);
		return globalIterate.add(sol);
	}
}
