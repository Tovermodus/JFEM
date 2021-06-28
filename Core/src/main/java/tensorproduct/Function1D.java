package tensorproduct;

public abstract class Function1D
{
	public abstract double value(double pos);

	public abstract double derivative(double pos);
	public static Function1D constFunction(double con)
	{
		return new Function1D()
		{
			@Override
			public double value(double pos)
			{
				return con;
			}

			@Override
			public double derivative(double pos)
			{
				return 0;
			}
		};
	}
}
