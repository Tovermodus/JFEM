package schwarz;

import linalg.Vector;
import linalg.VectorMultiplyable;
import org.jetbrains.annotations.NotNull;

public interface SubspaceCorrection<OT extends VectorMultiplyable>
{
	
	Vector apply(final AbstractSchwarz<?, ?, OT> schwarz, final @NotNull Vector globalIterate,
	             final @NotNull Vector globalRhs);
}
