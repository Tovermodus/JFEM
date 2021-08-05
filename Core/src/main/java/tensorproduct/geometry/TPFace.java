package tensorproduct.geometry;

import basic.*;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableSet;
import com.google.common.collect.Multimap;
import linalg.CoordinateComparator;
import linalg.CoordinateVector;
import org.jetbrains.annotations.NotNull;

import java.util.*;

public class TPFace implements FaceWithReferenceFace<TPCell, TPFace>,Comparable<TPFace>
{
	public final double otherCoordinate;
	final ImmutableList<Cell1D> cell1Ds;
	public final int flatDimension;
	
	Set<TPCell> cells;
	private final boolean isBoundaryFace;
	private final VectorFunction normal;
	
	TPFace(List<Cell1D> cell1Ds, int flatDimension, double otherCoordinate, boolean isBoundaryFace)
	{
		this.otherCoordinate = otherCoordinate;
		this.flatDimension = flatDimension;
		this.isBoundaryFace = isBoundaryFace;
		this.cell1Ds = ImmutableList.copyOf(cell1Ds);
		this.cells = new HashSet<>();
		this.normal = new VectorFunction()
		{
			@Override
			public int getRangeDimension()
			{
				return getDimension();
			}
			
			@Override
			public int getDomainDimension()
			{
				return getDimension();
			}
			
			@Override
			public CoordinateVector value(CoordinateVector pos)
			{
				CoordinateVector ret = new CoordinateVector(pos.getLength());
				ret.set(1,flatDimension);
				return ret;
			}
		};
	}
	public ImmutableList<Cell1D> getComponentCells()
	{
		return cell1Ds;
	}
	@Override
	public int getDimension()
	{
		return cell1Ds.size()+1;
	}
	
	@Override
	public ImmutableSet<TPCell> getCells()
	{
		return ImmutableSet.copyOf(cells);
	}
	
	@Override
	public boolean isBoundaryFace()
	{
		return isBoundaryFace;
	}
	
	@Override
	public VectorFunction getNormal()
	{
		return normal;
	}
	public boolean isNormalDownstream(CoordinateVector pos)
	{
		return getNormal().value(center()).inner(pos.sub(center()))>0;
	}
	public TPCell getUpStreamCell(CoordinateVector direction)
	{
		if(normal.value(center()).inner(direction)>0)
			return getNormalUpstreamCell();
		else
			return getNormalDownstreamCell();
	}
	public TPCell getDownStreamCell(CoordinateVector direction)
	{
		if(normal.value(center()).inner(direction)<0)
			return getNormalUpstreamCell();
		else
			return getNormalDownstreamCell();
	}
	
	@Override
	public TPCell getNormalDownstreamCell()
	{
		for(TPCell cell:cells)
		{
			if (cell.center().at(flatDimension) > otherCoordinate)
			{
				return cell;
			}
		}
		return null;
	}
	@Override
	public TPCell getNormalUpstreamCell()
	{
		for(TPCell cell:cells)
		{
			if (cell.center().at(flatDimension) < otherCoordinate)
			{
				return cell;
			}
		}
		return null;
	}
	
	@Override
	public CoordinateVector center()
	{
		CoordinateVector ret = new CoordinateVector(getDimension());
		int subd = 0;
		for(int d = 0; d < getDimension(); d++)
		{
			if(d == flatDimension)
				ret.set(otherCoordinate, d);
			else
				ret.set(cell1Ds.get(subd++).center(),d);
		}
		return ret;
	}
	
	@Override
	public boolean isOnFace(CoordinateVector pos)
	{
		int subd = 0;
		for(int d = 0; d < getDimension(); d++)
		{
			if(d == flatDimension)
			{
				if (!DoubleCompare.almostEqual(pos.at(flatDimension), otherCoordinate))
					return false;
			}
			else
			{
				if (!cell1Ds.get(subd++).isInCell(pos.at(d)))
					return false;
			}
		}
		return true;
		
	}
	
	@Override
	public List<TPFace> refine(Multimap<TPCell, TPCell> cellMap)
	{
		throw new UnsupportedOperationException();
	}

	@Override
	public String toString()
	{
		String ret = "";
		int subd = 0;
		for(int d = 0; d < getDimension(); d++)
		{
			if(d == flatDimension)
				ret = ret.concat(otherCoordinate+"");
			else
				ret = ret.concat("["+cell1Ds.get(subd).getStart()+", "+cell1Ds.get(subd++).getEnd()+
					"]");
			if(d<cell1Ds.size())
				ret = ret.concat("x");
		}
		if(isBoundaryFace)
			ret = ret.concat(" On Boundary");
		return ret;
	}
	@Override
	public int compareTo(@NotNull TPFace o)
	{
		if(o.getDimension() < getDimension())
			return -1;
		if(o.getDimension() > getDimension())
			return 1;
		if(o.flatDimension < flatDimension)
			return -1;
		if(o.flatDimension > flatDimension)
			return 1;
		if(o.otherCoordinate < otherCoordinate-PerformanceArguments.getInstance().doubleTolerance)
			return -1;
		if(o.otherCoordinate > otherCoordinate+PerformanceArguments.getInstance().doubleTolerance)
			return 1;
		if(o.isBoundaryFace && !isBoundaryFace)
			return 1;
		if(!o.isBoundaryFace && isBoundaryFace)
			return -1;
		if(o.getNormalDownstreamCell() != null && getNormalDownstreamCell() == null)
			return 1;
		if(o.getNormalDownstreamCell() != null && getNormalDownstreamCell() == null)
			return -1;
		if(o.getNormalUpstreamCell() != null && getNormalUpstreamCell() == null)
			return 1;
		if(o.getNormalUpstreamCell() != null && getNormalUpstreamCell() == null)
			return -1;
		return CoordinateComparator.comp(center().getEntries(), o.center().getEntries());
	}
	@Override
	public int hashCode()
	{
		int ret = 0;
		for(int i = 0; i < cell1Ds.size(); i++)
		{
			ret += Math.pow(7, 3*i) * DoubleCompare.doubleHash(cell1Ds.get(i).center());
			ret += Math.pow(7, 3*i+1) * DoubleCompare.doubleHash(cell1Ds.get(i).getStart());
			ret += Math.pow(7, 3*i+2) * DoubleCompare.doubleHash(cell1Ds.get(i).getEnd());
		}
		ret -= 141*DoubleCompare.doubleHash(otherCoordinate);
		ret*=(flatDimension+19);
		ret*=2;
		if(isBoundaryFace)
			ret+=1;
		if(getNormalDownstreamCell() != null)
			ret+= 1789821;
		if(getNormalUpstreamCell() != null)
			ret+= 31789821;
		return ret;
	}
	
	@Override
	public boolean equals(Object obj)
	{
		if (obj instanceof TPFace)
			return 0 == compareTo((TPFace) obj);
		else
			return false;
	}
	
	@Override
	public TPFace getReferenceFace()
	{
		List<Cell1D> cells1 = new ArrayList<>(getDimension()-1);
		for(int i = 0; i < getDimension()-1; i++)
		{
			cells1.add(new Cell1D(0, 1));
		}
		TPFace refFace =
			new TPFace(cells1, flatDimension, 0, isBoundaryFace);
		cells1 = new ArrayList<>(getDimension()-1);
		List<Cell1D> cells2 = new ArrayList<>(getDimension()-1);
		for (int i = 0; i < getDimension(); i++)
		{
			if(i != flatDimension)
			{
				cells1.add(new Cell1D(0, 1));
				cells2.add(new Cell1D(0, 1));
			}
			else
			{
				cells1.add(new Cell1D(0, 1));
				cells2.add(new Cell1D(-1, 0));
			}
		}
		TPCell downstreamCell = new TPCell(cells1);
		if(getNormalDownstreamCell() != null)
			refFace.cells.add(downstreamCell);
		TPCell upstreamCell = new TPCell(cells2);
		if(getNormalUpstreamCell() != null)
			refFace.cells.add(upstreamCell);
		return refFace;
	}
}