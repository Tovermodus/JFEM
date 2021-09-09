package tensorproduct.geometry;

import basic.DoubleCompare;
import basic.PerformanceArguments;
import org.jetbrains.annotations.NotNull;
import tensorproduct.QuadratureRule1D;

public class Cell1D implements Comparable<Cell1D>
{
	private final double start;
	private final double end;
	
	Cell1D(final double start, final double end)
	{
		this.start = start;
		this.end = end;
	}
	
	Cell1D(final Cell1D cell)
	{
		this(cell.start, cell.end);
	}
	
	public double[][] distributeQuadrature(final QuadratureRule1D quadratureRule)
	{
		final double[][] pointsWeights = new double[2][quadratureRule.length()];
		for (int i = 0; i < quadratureRule.length(); i++)
		{
			pointsWeights[0][i] = positionOnGrid(quadratureRule.getReferencePoints().get(i));
			pointsWeights[1][i] = quadratureRule.getReferenceWeights().get(i) * length();
		}
		return pointsWeights;
	}
	
	public double getStart()
	{
		return start;
	}
	
	public double getEnd()
	{
		return end;
	}
	
	public double length()
	{
		return (end - start);
	}
	
	public double center()
	{
		return (end + start) / 2;
	}
	
	public boolean isInCell(final double pos)
	{
		return (pos >= start - 100 * PerformanceArguments.getInstance().doubleTolerance &&
			        pos <= end + 100 * PerformanceArguments.getInstance().doubleTolerance);
	}
	
	public double positionOnReferenceCell(final double pos)
	{
		return (pos - start) / length();
	}
	
	public double positionOnGrid(final double pospp)
	{
		return pospp * length() + start;
	}
	
	@Override
	public String toString()
	{
		return "Cell: start" + start + ", end: " + end;
	}
	
	public void print()
	{
		System.out.println(this);
	}
	
	@Override
	public int compareTo(@NotNull final Cell1D o)
	{
		if (getStart() < o.getStart() - PerformanceArguments.getInstance().doubleTolerance)
			return -1;
		if (getStart() > o.getStart() + PerformanceArguments.getInstance().doubleTolerance)
			return 1;
		if (getEnd() < o.getEnd() - PerformanceArguments.getInstance().doubleTolerance)
			return -1;
		if (getEnd() > o.getEnd() + PerformanceArguments.getInstance().doubleTolerance)
			return 1;
		return 0;
	}
	
	@Override
	public int hashCode()
	{
		return 17 + DoubleCompare.doubleHash(start) + 39 * DoubleCompare.doubleHash(end);
	}
	
	@Override
	public boolean equals(final Object obj)
	{
		if (!(obj instanceof Cell1D))
			return false;
		return compareTo((Cell1D) obj) == 0;
	}
}
