package schwarz;

import basic.Cell;
import basic.Face;
import linalg.Matrix;

import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public abstract class UpFrontSchwarz<CT extends Cell<CT, FT>, FT extends Face<CT, FT>>
	implements MatrixSchwarz<CT, FT>
{
	List<Matrix> restrictionOperators;
	List<Matrix> localOperators;
	Matrix globalMatrix;
	
	public UpFrontSchwarz(final Matrix globalMatrix)
	{
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
			         .parallel()
			         .mapToObj(i -> restrictionOperators.get(i)
			                                            .mmMul(globalMatrix)
			                                            .mtMul(restrictionOperators.get(i)))
			         .collect(Collectors.toList());
	}
	
	public abstract Matrix buildRestrictionMatrix(int patch);
	
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
