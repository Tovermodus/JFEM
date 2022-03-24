package schwarz;

import linalg.DenseVector;
import linalg.Vector;
import linalg.VectorMultiplyable;

import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class AdditiveSubspaceCorrection<OT extends VectorMultiplyable>
	implements SubspaceCorrection<OT>
{
	private final double omega;
	
	public AdditiveSubspaceCorrection(final double omega)
	{
		this.omega = omega;
	}
	
	@Override
	public Vector solve(final AbstractSchwarz<?, ?, OT> schwarz, final Vector globalRhs)
	{
		final List<Vector> globalSolutionComponents
			= IntStream.range(0, schwarz.getPatchCount())
			           .parallel()
			           .mapToObj(i ->
			                     {
				                     final Vector localRhs
					                     = schwarz.getLocalVector(i,
					                                              globalRhs);
				                     final Vector localSol
					                     = schwarz.solveLocalSystem(i,
					                                                localRhs);
				                     final Vector globalSolComponent
					                     = schwarz.getGlobalVector(i,
					                                               localSol);
				                     return globalSolComponent;
			                     })
			           .collect(Collectors.toList());
		Vector ret = new DenseVector(globalRhs.mul(0));
		for (final Vector v : globalSolutionComponents)
			ret = ret.add(v.mul(omega));
		return ret;
	}
}
