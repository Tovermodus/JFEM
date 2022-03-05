package schwarz;

import basic.Cell;
import basic.Face;
import linalg.DenseMatrix;
import linalg.Matrix;
import linalg.Vector;

public interface MatrixSchwarz<CT extends Cell<CT, FT>, FT extends Face<CT, FT>>
	extends AbstractSchwarz<CT, FT, Matrix>
{
	@Override
	default Vector solveLocalSystem(final int patch, final Vector localVector)
	{
		final DenseMatrix localOp = new DenseMatrix(getLocalOperator(patch));
		return localOp.solve(localVector);
	}
}
