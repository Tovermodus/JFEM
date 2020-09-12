package tensorproduct;

import com.google.common.primitives.Doubles;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public enum QuadratureRule1D
{
	Gauss0(new double[]{}, //on [-1,1]
		new double[]{}),
	Gauss1(new double[]{2}, //on [-1,1]
		new double[]{0}),
	Gauss2(new double[]{1, 1}, //on [-1,1]
		new double[]{-Math.sqrt(1. / 3), Math.sqrt(1. / 3)}),
	Gauss3(new double[]{5. / 9, 8. / 9, 5. / 9}, //on [-1,1]
		new double[]{-Math.sqrt(3. / 5), 0, Math.sqrt(3. / 5)}),
	Gauss4(new double[]{
		(18. - Math.sqrt(30.)) / 36,
		(18. + Math.sqrt(30.)) / 36,
		(18. + Math.sqrt(30.)) / 36,
		(18. - Math.sqrt(30.)) / 36}, //on [-1,1]
		new double[]{
			-Math.sqrt(3. / 7 + 2. / 7 * Math.sqrt(6. / 5)),
			-Math.sqrt(3. / 7 - 2. / 7 * Math.sqrt(6. / 5)),
			Math.sqrt(3. / 7 - 2. / 7 * Math.sqrt(6. / 5)),
			Math.sqrt(3. / 7 + 2. / 7 * Math.sqrt(6. / 5))
		}),
	Gauss5(new double[]{
		(322.0 - 13 * Math.sqrt(70)) / 900,
		(322.0 + 13 * Math.sqrt(70)) / 900,
		128.0 / 225,
		(322.0 + 13 * Math.sqrt(70)) / 900,
		(322.0 - 13 * Math.sqrt(70)) / 900},
		new double[]{
			-Math.sqrt(5.0 + 2 * Math.sqrt(10.0 / 7)) / 3,
			-Math.sqrt(5.0 - 2 * Math.sqrt(10.0 / 7)) / 3,
			0,
			Math.sqrt(5.0 - 2 * Math.sqrt(10.0 / 7)) / 3,
			Math.sqrt(5.0 + 2 * Math.sqrt(10.0 / 7)) / 3}),
	GaussLobatto2(new double[]{1. / 2, 1. / 2}, new double[]{-1, 1}),
	GaussLobatto3(new double[]{1. / 3, 4. / 3, 1. / 3}, new double[]{-1, 0, 1}),
	GaussLobatto4(new double[]{1. / 6, 5. / 6, 5. / 6, 1. / 6},
		new double[]{-1, -Math.sqrt(1. / 5), Math.sqrt(1. / 5), 1}),
	GaussLobatto5(new double[]{1. / 10, 49. / 90, 32. / 45, 49. / 90, 1. / 10}, new double[]{
		-1, -Math.sqrt(3. / 7), 0, Math.sqrt(3. / 7), 1
	}),
	GaussLobatto6(new double[]{
		1. / 15,
		(14. - Math.sqrt(7)) / 30,
		(14. + Math.sqrt(7)) / 30,
		(14. + Math.sqrt(7)) / 30,
		(14. - Math.sqrt(7)) / 30,
		1. / 15
	}, new double[]{
		-1,
		-Math.sqrt(1. / 3 + 2 * Math.sqrt(7) / 21),
		-Math.sqrt(1. / 3 - 2 * Math.sqrt(7) / 21),
		+Math.sqrt(1. / 3 - 2 * Math.sqrt(7) / 21),
		+Math.sqrt(1. / 3 + 2 * Math.sqrt(7) / 21),
		1
	}
	);
	private final double[] referenceWeights;
	private final double[] referencePoints;
	
	public final int length()
	{
		return referencePoints.length;
	}
	
	public final List<Double> getReferenceWeights()
	{
		return Collections.unmodifiableList(Doubles.asList(referenceWeights));
	}
	
	public final List<Double> getReferencePoints()
	{
		return Collections.unmodifiableList(Doubles.asList(referencePoints));
	}
	
	QuadratureRule1D(double[] referenceWeights, double[] referencePoints)
	{
		this.referenceWeights = referenceWeights;
		this.referencePoints = referencePoints;
		for (int i = 0; i < referencePoints.length; i++)
		{
			this.referencePoints[i] = referencePoints[i] / 2 + 0.5;
			this.referenceWeights[i] = referenceWeights[i] / 2;
		}
	}
	
}
