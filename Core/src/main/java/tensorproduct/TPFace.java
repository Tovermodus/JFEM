package tensorproduct;

import basic.Face;
import basic.FaceWithReferenceFace;
import basic.VectorFunction;
import com.google.common.collect.Multimap;
import linalg.CoordinateComparator;
import linalg.CoordinateVector;
import org.jetbrains.annotations.NotNull;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

public class TPFace implements FaceWithReferenceFace<TPCell, TPFace, TPEdge>,Comparable<TPFace>
{
	double otherCoordinate;
	List<Cell1D> cell1Ds;
	int flatDimension;
	
	private final Set<TPCell> cells;
	private final Set<TPEdge> edges;
	private boolean isBoundaryFace;
	private final VectorFunction normal;
	
	public TPFace(List<Cell1D> cell1Ds, int flatDimension, double otherCoordinate, boolean isBoundaryFace)
	{
		this.otherCoordinate = otherCoordinate;
		this.flatDimension = flatDimension;
		this.isBoundaryFace = isBoundaryFace;
		this.cell1Ds = new ArrayList<>();
		for(Cell1D c: cell1Ds)
			this.cell1Ds.add(new Cell1D(c));
		this.edges = new TreeSet<>();
		this.cells = new TreeSet<>();
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
	
	@Override
	public int getDimension()
	{
		return cell1Ds.size()+1;
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
	
	
	public List<Cell1D> getCell1Ds()
	{
		return cell1Ds;
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
	@Override
	public int hashCode()
	{
		int ret = 0;
		for(int i = 0; i < cell1Ds.size(); i++)
		{
			ret += Math.pow(7, 3*i) * cell1Ds.get(i).center();
			ret += Math.pow(7, 3*i+1) * cell1Ds.get(i).getStart();
			ret += Math.pow(7, 3*i+2) * cell1Ds.get(i).getEnd();
		}
		ret -= 141*otherCoordinate;
		ret*=(flatDimension+19);
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
	public void addEdge(TPEdge tpEdge)
	{
		if(edges.add(tpEdge))
			tpEdge.addFace(this);
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
		refFace.addCell(downstreamCell);
		TPCell upstreamCell = new TPCell(cells2);
		refFace.addCell(upstreamCell);
		return refFace;
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
