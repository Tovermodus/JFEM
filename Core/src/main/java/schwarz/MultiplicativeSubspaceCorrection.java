package schwarz;

import linalg.DenseVector;
import linalg.Vector;
import linalg.VectorMultiplyable;
import org.jetbrains.annotations.NotNull;

public class MultiplicativeSubspaceCorrection<OT extends VectorMultiplyable>
	implements SubspaceCorrection<OT>
{
	private final double omega;
	
	public MultiplicativeSubspaceCorrection(final double omega)
	{
		this.omega = omega;
	}
	
	@Override
	public Vector solve(final AbstractSchwarz<?, ?, OT> schwarz, final Vector globalRhs)
	{
		throw new UnsupportedOperationException("not implemented yet");
	}
	
	@Override
	public Vector apply(final AbstractSchwarz<?, ?, OT> schwarz,
	                    @NotNull final Vector globalIterate,
	                    @NotNull final Vector globalRhs)
	{
		Vector iterate = new DenseVector(globalIterate);
		for (int i = 0; i < schwarz.getPatchCount(); i++)
		{
			final Vector globalResidual = globalRhs.sub(schwarz.getGlobalOperator()
			                                                   .mvMul(iterate));
			//System.out.println("SUB " + globalResidual.euclidianNorm());
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
