package tensorproduct;

public abstract class Function1D
{
	public abstract double value(double pos);

	public abstract double derivative(double pos);
	public static Function1D constFunction(double con)
	{
		Function1D constF = new Function1D()
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
		return constF;
	}
}
