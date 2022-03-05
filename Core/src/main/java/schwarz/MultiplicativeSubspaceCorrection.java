package schwarz;

import linalg.DenseVector;
import linalg.Vector;
import linalg.VectorMultiplyable;
import org.jetbrains.annotations.NotNull;

public class MultiplicativeSubspaceCorrection<OT extends VectorMultiplyable>
	implements SubspaceCorrection<OT>
{
	private final double omega;
	
	public MultiplicativeSubspaceCorrection(double omega)
	{
		this.omega = omega;
	}
	@Override
	public Vector apply(final AbstractSchwarz<?, ?, OT> schwarz,
	                    @NotNull final Vector globalIterate,
	                    @NotNull final Vector globalRhs)
	{
		Vector iterate = new DenseVector(globalIterate);
		for (int i = 0; i < schwarz.getPatchCount(); i++)
		{
			final Vector globalResidual = globalRhs.sub(schwarz.getGlobalOperator().mvMul(globalIterate));
			final Vector localRes = schwarz.getLocalVector(i,
			                                       globalResidual);
			final Vector localSol = schwarz.solveLocalSystem(i,
			                                         localRes);
			final Vector globalSolComponent
				= schwarz.getGlobalVector(i,
				                  localSol);
			iterate = iterate.add(globalSolComponent.mul(omega));
		}
		return iterate;
	}
}
