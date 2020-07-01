package tensorproduct;

import com.google.common.primitives.Doubles;

import java.util.Arrays;
import java.util.Collections;
import java.util.List;

public enum QuadratureRule1D
{
	Gauss3(new double[]{5./9,8./9,5./9}, //on [-1,1]
		new double[]{-Math.sqrt(3./5),0,-Math.sqrt(3./5)}),
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
			Math.sqrt(5.0 + 2 * Math.sqrt(10.0 / 7)) / 3});
	private final double [] referenceWeights;
	private final double [] referencePoints;
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
	
	QuadratureRule1D(double [] referenceWeights, double [] referencePoints)
	{
		this.referenceWeights = referenceWeights;
		this.referencePoints = referencePoints;
		for(int i = 0; i < referencePoints.length; i++)
		{
			this.referencePoints[i] = referencePoints[i]/2+0.5;
			this.referenceWeights[i] = referenceWeights[i]/2;
		}
	}
	
}
