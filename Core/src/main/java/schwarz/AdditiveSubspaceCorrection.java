package schwarz;

import linalg.DenseVector;
import linalg.Vector;
import linalg.VectorMultiplyable;
import org.jetbrains.annotations.NotNull;

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
	public Vector apply(final AbstractSchwarz<?, ?, OT> schwarz,
	                    @NotNull final Vector globalIterate,
	                    @NotNull final Vector globalRhs)
	{
		final Vector globalResidual = globalRhs.sub(schwarz.getGlobalOperator()
		                                                   .mvMul(globalIterate));
		final List<Vector> globalSolutionComponents = IntStream.range(0, schwarz.getPatchCount())
		                                                       .parallel()
		                                                       .mapToObj(i ->
		                                                                 {
			                                                                 final Vector localRes
				                                                                 = schwarz.getLocalVector(
				                                                                 i,
				                                                                 globalResidual);
			                                                                 final Vector localSol
				                                                                 = schwarz.solveLocalSystem(
				                                                                 i,
				                                                                 localRes);
			                                                                 final Vector globalSolComponent
				                                                                 = schwarz.getGlobalVector(
				                                                                 i,
				                                                                 localSol);
			                                                                 return globalSolComponent;
		                                                                 })
		                                                       .collect(Collectors.toList());
		Vector ret = new DenseVector(globalIterate);
		for (final Vector v : globalSolutionComponents)
			ret = ret.add(v.mul(omega));
		return ret;
	}
}
