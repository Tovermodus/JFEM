package schwarz;

import linalg.DenseVector;
import linalg.Vector;
import linalg.VectorMultiplyable;
import org.jetbrains.annotations.NotNull;

public class MultiplicativeSubspaceCorrection<OT extends VectorMultiplyable>
	implements SubspaceCorrection<OT>
{
	@Override
	public Vector solve(final AbstractSchwarz<?, ?, OT> schwarz, final Vector globalRhs)
	{
		Vector iterate = schwarz.getGlobalVector(0, schwarz.solveLocalSystem(0,
		                                                                     schwarz.getLocalVector(0,
		                                                                                            globalRhs)));
		for (int i = 1; i < schwarz.getPatchCount(); i++)
		{
			
			final Vector globalResidual = globalRhs.sub(schwarz.getGlobalOperator()
			                                                   .mvMul(iterate));
			final Vector localRes = schwarz.getLocalVector(i,
			                                               globalResidual);
			final Vector localSol = schwarz.solveLocalSystem(i,
			                                                 localRes);
			final Vector globalSolComponent
				= schwarz.getGlobalVector(i,
				                          localSol);
			iterate = iterate.add(globalSolComponent);
		}
		return iterate;
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
			iterate = iterate.add(globalSolComponent);
		}
		return iterate;
	}
}
