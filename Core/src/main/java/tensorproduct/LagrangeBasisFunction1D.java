package tensorproduct;

import basic.DoubleCompare;
import tensorproduct.geometry.Cell1D;

public class LagrangeBasisFunction1D
	extends Function1D
{
	private final Cell1D cell;
	private final int polynomialDegree;
	private final int localFunctionNumber;
	private final double degreeOfFreedom;
	
	public LagrangeBasisFunction1D(final int polynomialDegree, final int localFunctionNumber, final Cell1D cell)
	{
		this.cell = cell;
		this.localFunctionNumber = localFunctionNumber;
		this.polynomialDegree = polynomialDegree;
		this.degreeOfFreedom = cell.positionOnGrid(1.0 / polynomialDegree * localFunctionNumber); //equidistant
	}
	
	public LagrangeBasisFunction1D(final int polynomialDegree, final double degreeOfFreedom, final Cell1D cell)
	{
		this.cell = cell;
		if (!cell.isInCell(degreeOfFreedom))
			throw new IllegalArgumentException("degree of freedom is not in given cell");
		if (!DoubleCompare.isInteger(cell.positionOnReferenceCell(degreeOfFreedom) * polynomialDegree))
			throw new IllegalArgumentException("Identification not possible. This should be int" + cell.positionOnReferenceCell(
				degreeOfFreedom) * polynomialDegree);
		localFunctionNumber = getLocalFunctionNumberOnReferenceCell(polynomialDegree,
		                                                            cell.positionOnReferenceCell(degreeOfFreedom));
		this.polynomialDegree = polynomialDegree;
		this.degreeOfFreedom = degreeOfFreedom; //equidistant
	}
	
	protected static int getLocalFunctionNumberOnReferenceCell(final int polynomialDegree,
	                                                           final double degreeOfFreedom)
	{
		final int localFunctionNumber;
		localFunctionNumber =
			(int) Math.round(degreeOfFreedom * polynomialDegree);
		return localFunctionNumber;
	}
	
	public int getLocalFunctionNumber()
	{
		return localFunctionNumber;
	}
	
	public int getPolynomialDegree()
	{
		return polynomialDegree;
	}
	
	public static double valueOnReferenceCell(final double pos,
	                                          final int localFunctionNumber,
	                                          final int polynomialDegree)
	{
		switch (polynomialDegree)
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
				for (int i = 0; i <= polynomialDegree; i++)
				{
					if (i != localFunctionNumber)
					{
						ret
							*= (pos - 1. * i / polynomialDegree) / (1. * localFunctionNumber / polynomialDegree - 1. * i / polynomialDegree);
					}
				}
				return ret;
		}
		return 0;
	}
	
	public double valueOnReferenceCell(final double pos)
	{
		return valueOnReferenceCell(pos, this.localFunctionNumber, this.polynomialDegree);
	}
	
	public Cell1D getCell()
	{
		return cell;
	}
	
	public static double derivativeOnReferenceCell(final double pos,
	                                               final int localFunctionNumber,
	                                               final int polynomialDegree)
	{
		switch (polynomialDegree)
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
				for (int j = 0; j <= polynomialDegree; j++)
				{
					if (j != localFunctionNumber)
					{
						double ret =
							1. / (1. * localFunctionNumber / polynomialDegree - 1. * j / polynomialDegree);
						for (int i = 0; i <= polynomialDegree; i++)
						{
							if (i != localFunctionNumber && i != j)
							{
								ret
									*= (pos - 1. * i / polynomialDegree) / (1. * localFunctionNumber / polynomialDegree - 1. * i / polynomialDegree);
							}
						}
						derret += ret;
					}
				}
				return derret;
		}
		return 0;
	}
	
	public double derivativeOnReferenceCell(final double pos)
	{
		return derivativeOnReferenceCell(pos, localFunctionNumber, polynomialDegree);
	}
	
	@Override
	public double value(final double pos)
	{
		if (cell.isInCell(pos))
			return valueOnReferenceCell(cell.positionOnReferenceCell(pos));
		else
			return 0;
	}
	
	@Override
	public double derivative(final double pos)
	{
		if (cell.isInCell(pos))
			return derivativeOnReferenceCell(cell.positionOnReferenceCell(pos)) / cell.length();
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
