package basic;

import linalg.DoubleTensor;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;

public abstract class ScalarFunction
{
	public abstract double value(DoubleTensor pos);

	public abstract DoubleTensor derivative(DoubleTensor pos);
	public static ScalarFunction oneFunction()
	{
		return constantFunction(1);
	}
	public static ScalarFunction constantFunction(double constant)
	{
		ScalarFunction oneF = new ScalarFunction()
		{
			@Override
			public double value(DoubleTensor pos)
			{
				return constant;
			}

			@Override
			public DoubleTensor derivative(DoubleTensor pos)
			{
				return new DoubleTensor((int)pos.size());
			}
		};
		return oneF;
	}
	public void plot(int pointres, String filename)
	{
		double[][] values = new double[pointres][pointres];
		for (int k = 0; k < pointres; k++)
		{
			for (int j = 0; j < pointres; j++)
			{
				DoubleTensor pos = DoubleTensor.vectorFromValues(1.0 / (pointres - 1) * k,
					1.0 / (pointres - 1) * j);

				values[k][j] = value(pos);
			}
		}
		try
		{
			BufferedWriter plotWriter = new BufferedWriter(new FileWriter(filename));
			for (int i = 0; i < pointres; i++)
			{
				for (int j = 0; j < pointres; j++)
				{
					plotWriter.write(Double.toString(values[i][j]) + " ");
				}
				plotWriter.newLine();
			}
			plotWriter.flush();
			plotWriter.close();
		} catch (
			IOException e)
		{
			e.printStackTrace();
		}
	}

}
