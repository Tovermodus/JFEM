package tensorproduct;

import basic.CellIntegral;
import basic.Function;
import basic.ScalarShapeFunction;
import linalg.CoordinateVector;
import linalg.Vector;

import java.util.List;
import java.util.function.ToDoubleFunction;

public class TPEdgeIntegral
{
	public static double integrateNonTensorProduct(ToDoubleFunction<CoordinateVector> eval, TPEdge e)
	{
		double v = TPCellIntegral.integrateNonTensorProduct(x ->
		{
			double[] point = new double[3];
			int subd = 0;
			for (int j = 0; j < 3; j++)
			{
				if (j == e.tangentialDimension)
					point[j] = x.at(e.tangentialDimension);
				else
					point[j] = e.otherCoordinates[subd++];
			}
			return eval.applyAsDouble(CoordinateVector.fromValues(point));
		}, List.of(e.cell));
		return v;
	}
}