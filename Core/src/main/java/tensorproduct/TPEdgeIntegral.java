package tensorproduct;

import linalg.CoordinateVector;

import java.util.List;
import java.util.function.ToDoubleFunction;

public class TPEdgeIntegral
{
	public static double integrateNonTensorProduct(ToDoubleFunction<CoordinateVector> eval, TPEdge e,
	                                               QuadratureRule1D quadratureRule)
	{
		return TPCellIntegral.integrateNonTensorProduct(x ->
		{
			double[] point = new double[3];
			int subd = 0;
			for (int j = 0; j < 3; j++)
			{
				if (j == e.tangentialDimension)
					point[j] = x.at(0);
				else
					point[j] = e.otherCoordinates[subd++];
			}
			return eval.applyAsDouble(CoordinateVector.fromValues(point));
		}, List.of(e.cell), quadratureRule);
	}
}