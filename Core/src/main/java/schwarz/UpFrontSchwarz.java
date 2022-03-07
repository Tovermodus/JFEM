package schwarz;

import basic.Cell;
import basic.Face;
import linalg.Matrix;

import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public abstract class UpFrontSchwarz<CT extends Cell<CT, FT>, FT extends Face<CT, FT>>
	extends AbstractSchwarz<CT, FT, Matrix>
{
	List<RestrictionMatrix> restrictionOperators;
	List<Matrix> localOperators;
	Matrix globalMatrix;
	
	public UpFrontSchwarz(final Matrix globalMatrix, final SubspaceCorrection<Matrix> subspaceCorrection,
	                      final SystemSolver<Matrix> solver)
	{
		super(subspaceCorrection, solver);
		this.globalMatrix = globalMatrix;
	}
	
	public void build()
	{
		restrictionOperators =
			IntStream.range(0, getPatchCount())
			         .parallel()
			         .mapToObj(this::buildRestrictionMatrix)
			         .collect(
				         Collectors.toList());
		localOperators =
			IntStream.range(0, getPatchCount())
			         .mapToObj(i ->
				                   restrictionOperators.get(i)
				                                       .selectFrom(globalMatrix))
			         .collect(Collectors.toList());
	}
	
	public abstract RestrictionMatrix buildRestrictionMatrix(int patch);
	
	@Override
	public Matrix getGlobalOperator()
	{
		return globalMatrix;
	}
	
	@Override
	public Matrix getRestrictionOperator(final int patch)
	{
		return restrictionOperators.get(patch);
	}
	
	@Override
	public Matrix getLocalOperator(final int patch)
	{
		return localOperators.get(patch);
	}
}
