package tensorproduct;

import basic.Face;
import basic.ShapeFunction;
import basic.VectorFunction;
import com.google.common.collect.Multimap;
import basic.Cell;
import linalg.CoordinateComparator;
import linalg.CoordinateVector;
import linalg.Vector;
import org.jetbrains.annotations.NotNull;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

public class TPFace implements Face<TPCell, TPFace>,Comparable<TPFace>
{
	double otherCoordinate;
	List<Cell1D> cell1Ds;
	int flatDimension;
	private Set<TPCell> cells;
	private boolean isBoundaryFace;
	private VectorFunction normal;
	
	public TPFace(List<Cell1D> cell1Ds, int flatDimension, double otherCoordinate, boolean isBoundaryFace)
	{
		this.otherCoordinate = otherCoordinate;
		this.flatDimension = flatDimension;
		this.isBoundaryFace = isBoundaryFace;
		this.cell1Ds = cell1Ds;
		this.cells = new TreeSet<>();
		this.normal = new VectorFunction()
		{
			@Override
			public int getDomainDimension()
			{
				return getDimension()+1;
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
	
	@Override
	public int getDimension()
	{
		return cell1Ds.size();
	}
	
	@Override
	public Set<TPCell> getCells()
	{
		return cells;
	}
	
	@Override
	public void setBoundaryFace(boolean boundaryFace)
	{
		isBoundaryFace = boundaryFace;
	}
	
	@Override
	public boolean isBoundaryFace()
	{
		return isBoundaryFace;
	}
	
	@Override
	public void addCell(TPCell cell)
	{
		if(cells.add(cell))
			cell.addFace(this);
	}
	
	
	@Override
	public VectorFunction getNormal()
	{
		return normal;
	}
	public TPCell getUpStreamCell(CoordinateVector pos, CoordinateVector direction)
	{
		if(normal.value(pos).inner(direction)>0)
			return getNormalUpstreamCell(pos);
		else
			return getNormalDownstreamCell(pos);
	}
	public TPCell getDownStreamCell(CoordinateVector pos, CoordinateVector direction)
	{
		if(normal.value(pos).inner(direction)<0)
			return getNormalUpstreamCell(pos);
		else
			return getNormalDownstreamCell(pos);
	}
	
	@Override
	public TPCell getNormalDownstreamCell(CoordinateVector pos)
	{
		for(TPCell cell:cells)
		{
			if(cell.isInCell(pos))
			{
				if (cell.center().at(flatDimension) > otherCoordinate)
				{
					return cell;
				}
			}
		}
		return null;
	}
	@Override
	public TPCell getNormalUpstreamCell(CoordinateVector pos)
	{
		for(TPCell cell:cells)
		{
			if(cell.isInCell(pos))
			{
				if (cell.center().at(flatDimension) < otherCoordinate)
				{
					return cell;
				}
			}
		}
		return null;
	}
	
	@Override
	public CoordinateVector center()
	{
		CoordinateVector ret = new CoordinateVector(getDimension()+1);
		int subd = 0;
		for(int d = 0; d < getDimension()+1; d++)
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
		for(int d = 0; d < getDimension()+1; d++)
		{
			if(d == flatDimension)
			{
				if (pos.at(flatDimension) != otherCoordinate)
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

	public int getFlatDimension()
	{
		return flatDimension;
	}

	public double getOtherCoordinate()
	{
		return otherCoordinate;
	}
	
	@Override
	public String toString()
	{
		String ret = "";
		int subd = 0;
		for(int d = 0; d < getDimension()+1; d++)
		{
			if(d == flatDimension)
				ret = ret.concat(otherCoordinate+"");
			else
				ret = ret.concat("["+cell1Ds.get(subd).getStart()+", "+cell1Ds.get(subd++).getEnd()+
					"]");
			if(d<cell1Ds.size())
				ret = ret.concat("x");
		}
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
		if(o.otherCoordinate < otherCoordinate)
			return -1;
		if(o.otherCoordinate > otherCoordinate)
			return 1;
		return CoordinateComparator.comp(center().getEntries(), o.center().getEntries());
	}
//	@Override
//	public List<Face<TPCell, TPShapeFunction>> refine(Multimap<TPCell, TPCell> cellMap)
//	{
//		assert(cells.size()<2);
//		List<Face<TPCell, TPShapeFunction>> refinedFaces = new ArrayList<>();
//		TPFace face1 = new TPFace(new Cell1D(cell1d.getStart(),cell1d.center()),otherCoordinate,
//		normaldirection);
//		TPFace face2 = new TPFace(new Cell1D(cell1d.center(), cell1d.getEnd()),otherCoordinate,normaldirection);
//		for(TPCell cell:cells)
//		{
//			for (TPCell  refinedCell : cellMap.get(cell))
//			{
//				if(refinedCell.isInCell(face1.center()))
//				{
//					refinedCell.getFaces().add(face1);
//					face1.cells.add(refinedCell);
//				}
//				if(refinedCell.isInCell(face2.center()))
//				{
//					refinedCell.getFaces().add(face2);
//					face2.cells.add(refinedCell);
//				}
//			}
//		}
//		face1.setBoundaryFace(isBoundaryFace());
//		face2.setBoundaryFace(isBoundaryFace());
//		refinedFaces.add(face1);
//		refinedFaces.add(face2);
//		return refinedFaces;
//
//	}
//
//



}
