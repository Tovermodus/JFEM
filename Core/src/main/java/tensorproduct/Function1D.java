package tensorproduct;

public abstract class Function1D
{
	public abstract double value(double pos);

	public abstract double derivative(double pos);
	public static Function1D oneFunction()
	{
		Function1D oneF = new Function1D()
		{
			@Override
			public double value(double pos)
			{
				return 1;
			}

			@Override
			public double derivative(double pos)
			{
				return 0;
			}
		};
		return oneF;
	}
}
