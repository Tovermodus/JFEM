package tensorproduct;

import org.jetbrains.annotations.NotNull;

public class LagrangeBasisFunction1D extends Function1D
{
	private final Cell1D cell;
	private final int polynomialDegree;
	private final int localFunctionNumber;
	private final double degreeOfFreedom;
	
	public LagrangeBasisFunction1D(int polynomialDegree, int localFunctionNumber, Cell1D cell)
	{
		this.cell = cell;
		this.localFunctionNumber = localFunctionNumber;
		this.polynomialDegree = polynomialDegree;
		this.degreeOfFreedom = cell.positionOnGrid(1.0 / polynomialDegree * localFunctionNumber); //equidistant
	}
	
	public LagrangeBasisFunction1D(int polynomialDegree, double degreeOfFreedom, Cell1D cell)
	{
		this.cell = cell;
		if (!cell.isInCell(degreeOfFreedom))
			throw new IllegalArgumentException("degree of freedom is not in given cell");
		if (Math.abs(cell.positionOnReferenceCell(degreeOfFreedom) * polynomialDegree - (int) (cell.positionOnReferenceCell(degreeOfFreedom) * polynomialDegree + 1e-10)) > 1e-13)
			throw new IllegalArgumentException("Identification not possible " + cell.positionOnReferenceCell(degreeOfFreedom) * polynomialDegree + " " + (int) (cell.positionOnReferenceCell(degreeOfFreedom) * polynomialDegree + 1e-10));
		this.localFunctionNumber =
			(int) (cell.positionOnReferenceCell(degreeOfFreedom) * polynomialDegree + 1e-12);
		this.polynomialDegree = polynomialDegree;
		this.degreeOfFreedom = degreeOfFreedom; //equidistant
	}
	public int getLocalFunctionNumber()
	{
		return localFunctionNumber;
	}
	public int getPolynomialDegree()
	{
		return polynomialDegree;
	}
	public double valueOnReferenceCell(double pos)
	{
		switch (this.polynomialDegree)
		{
			case 1:
				switch (localFunctionNumber)
				{
					case 0:
						return 1.0 - pos;
					case 1:
						return pos;
				}
				break;
			case 2:
				switch (localFunctionNumber)
				{
					case 0:
						return 1.0 - 3.0 * pos + 2.0 * pos * pos;
					case 1:
						return (pos - pos * pos) * 4.0;
					case 2:
						return 2.0 * pos * pos - pos;
				}
				break;
			default:
				double ret = 1;
				for (int i = 0; i <= this.polynomialDegree; i++)
				{
					if (i != localFunctionNumber)
					{
						ret *= (pos - 1. * i / this.polynomialDegree) / (1. * localFunctionNumber / this.polynomialDegree - 1. * i / this.polynomialDegree);
					}
				}
				return ret;
		}
		return 0;
	}
	
	public double derivativeOnReferenceCell(double pos)
	{
		
		switch (this.polynomialDegree)
		{
			case 1:
				switch (localFunctionNumber)
				{
					case 0:
						return -1.0;
					case 1:
						return 1.0;
				}
				break;
			case 2:
				switch (localFunctionNumber)
				{
					case 0:
						return -3.0 + 4.0 * pos;
					case 1:
						return (1.0 - 2.0 * pos) * 4.0;
					case 2:
						return 4.0 * pos - 1.0;
				}
				break;
			default:
				double derret = 0;
				for (int j = 0; j <= this.polynomialDegree; j++)
				{
					if (j != localFunctionNumber)
					{
						double ret =
							1. / (1. * localFunctionNumber / this.polynomialDegree - 1. * j / this.polynomialDegree);
						for (int i = 0; i <= this.polynomialDegree; i++)
						{
							if (i != localFunctionNumber && i != j)
							{
								ret *= (pos - 1. * i / this.polynomialDegree) / (1. * localFunctionNumber / this.polynomialDegree - 1. * i / this.polynomialDegree);
							}
						}
						derret += ret;
					}
				}
				return derret;
		}
		return 0;
	}
	
	@Override
	public double value(double pos)
	{
		if (cell.isInCell(pos))
			return valueOnReferenceCell(cell.positionOnReferenceCell(pos));
		else
			return 0;
	}
	
	@Override
	public double derivative(double pos)
	{
		if (cell.isInCell(pos))
			return derivativeOnReferenceCell(cell.positionOnReferenceCell(pos)) * cell.jacobiDeterminant(pos);
		return 0;
	}
	
	public void print()
	{
		System.out.println("Shapefunction: Degree " + polynomialDegree + ", Local Number " + localFunctionNumber +
			", on Cell:\n\t\t");
		cell.print();
	}
	
	public double getDegreeOfFreedom()
	{
		return degreeOfFreedom;
	}
	
}
