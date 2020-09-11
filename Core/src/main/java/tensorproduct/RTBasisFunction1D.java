package tensorproduct;

public class RTBasisFunction1D extends Function1D
{
	private final Cell1D cell;
	private final int polynomialDegree;
	private int localFunctionNumber;
	private double degreeOfFreedom;
	QuadratureRule1D quad;
	
	public RTBasisFunction1D(int polynomialDegree, int localFunctionNumber, Cell1D cell, boolean higherDegree)
	{
		this.cell = cell;
		this.localFunctionNumber = localFunctionNumber;
		this.polynomialDegree = polynomialDegree;
		this.quad = obtainquadPoints(polynomialDegree, higherDegree);
		this.degreeOfFreedom = cell.positionOnGrid(quad.getReferencePoints().get(localFunctionNumber));
	}
	
	public RTBasisFunction1D(int polynomialDegree, double degreeOfFreedom, Cell1D cell, boolean higherDegree)
	{
		this.cell = cell;
		if (!cell.isInCell(degreeOfFreedom))
			throw new IllegalArgumentException("degree of freedom is not in given cell");
		this.quad = obtainquadPoints(polynomialDegree, higherDegree);
		boolean pointFound = false;
		for (int i = 0; i < quad.length(); i++)
		{
			if(Math.abs(degreeOfFreedom - cell.positionOnGrid(quad.getReferencePoints().get(i)))< 1e-12)
			{
				this.localFunctionNumber = i;
				this.degreeOfFreedom = quad.getReferencePoints().get(i);
				pointFound = true;
				break;
			}
		}
		if (!pointFound)
			throw new IllegalArgumentException("Identification not possible ");
		this.polynomialDegree = polynomialDegree;
	}
	private QuadratureRule1D obtainquadPoints(int polynomialDegree, boolean higherDegree)
	{
		
		if(higherDegree)
			switch (polynomialDegree)
			{
				case 1:
					return QuadratureRule1D.GaussLobatto2;
				case 2:
					return QuadratureRule1D.GaussLobatto3;
				case 3:
					return QuadratureRule1D.GaussLobatto4;
				case 4:
					return QuadratureRule1D.GaussLobatto5;
				case 5:
					return QuadratureRule1D.GaussLobatto6;
				default:
					throw new IllegalArgumentException("bad polynomial degree");
			}
		switch (polynomialDegree)
		{
			case 1:
				return QuadratureRule1D.Gauss1;
			case 2:
				return QuadratureRule1D.Gauss2;
			case 3:
				return QuadratureRule1D.Gauss3;
			case 4:
				return QuadratureRule1D.Gauss4;
			case 5:
				return QuadratureRule1D.Gauss5;
			default:
				throw new IllegalArgumentException("bad polynomial degree");
		}
	}
	public double valueOnReferenceCell(double pos)
	{
		if(!cell.isInCell(pos))
			return 0;
		double ret = 1;
		for (int i = 0; i < quad.length(); i++)
		{
			if (i != localFunctionNumber)
			{
				ret *= (pos - quad.getReferencePoints().get(i)) / (quad.getReferencePoints().get(localFunctionNumber) -  quad.getReferencePoints().get(i));
			}
		}
		return ret;
		
	}
	
	public double derivativeOnReferenceCell(double pos)
	{
		if (!cell.isInCell(pos))
			return 0;
		double derret = 0;
		for (int j = 0; j < quad.length(); j++)
		{
			if (j != localFunctionNumber)
			{
				double ret =
					1. / (quad.getReferencePoints().get(localFunctionNumber) - quad.getReferencePoints().get(j));
				for (int i = 0; i < quad.length(); i++)
				{
					if (i != localFunctionNumber && i != j)
					{
						ret *= (pos - quad.getReferencePoints().get(i)) /
							(quad.getReferencePoints().get(localFunctionNumber) - quad.getReferencePoints().get(i));
					}
				}
				derret += ret;
			}
		}
		return derret;
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
