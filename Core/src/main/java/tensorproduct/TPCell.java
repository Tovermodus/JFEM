package tensorproduct;

import basic.Cell;
import basic.CellWithReferenceCell;
import basic.ShapeFunction;
import basic.VectorFunction;
import linalg.CoordinateComparator;
import linalg.CoordinateMatrix;
import linalg.CoordinateVector;
import org.jetbrains.annotations.NotNull;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

public class TPCell implements CellWithReferenceCell<TPCell, TPFace, TPEdge>
{
	List<Cell1D> cell1Ds;
	Set<TPFace> faces;
	Set<TPEdge> edges;
	boolean refined;
	
	public TPCell(List<Cell1D> cell1Ds)
	{
		this.cell1Ds = new ArrayList<>();
		for(Cell1D c: cell1Ds)
			this.cell1Ds.add(new Cell1D(c));
		this.faces = new TreeSet<>();
		this.edges = new TreeSet<>();
		this.refined = false;
	}
	@Override
	public int getDimension()
	{
		return cell1Ds.size();
	}
	
	public List<Cell1D> getCell1Ds()
	{
		return cell1Ds;
	}
	
	@Override
	public Set<TPFace> getFaces()
	{
		return faces;
	}
	
	@Override
	public boolean isRefined()
	{
		return refined;
	}
	
	@Override
	public void setRefined(boolean refined)
	{
		refined = refined;
	}
	@Override
	public void addFace(TPFace face)
	{
		if(faces.add(face))
			face.addCell(this);
	}
	
	@Override
	public boolean isInCell(CoordinateVector pos)
	{
		for(int d = 0; d < cell1Ds.size(); d++)
		{
			if (!cell1Ds.get(d).isInCell(pos.at(d)))
				return false;
		}
		return true;
	}
	
	@Override
	public CoordinateVector center()
	{
		CoordinateVector ret = new CoordinateVector(getDimension());
		for(int i = 0; i < ret.getLength(); i++)
			ret.set(cell1Ds.get(i).center(),i);
		return ret;
	}
	
	@Override
	public VectorFunction getOuterNormal(TPFace face)
	{
		if(!faces.contains(face))
			throw new IllegalArgumentException("face does not belong to cell");
		boolean invertNormal = (center().sub(face.center())).inner(face.getNormal().value(face.center()))>0;
		if(invertNormal)
			return new VectorFunction()
			{
				@Override
				public int getRangeDimension()
				{
					return face.getNormal().getRangeDimension();
				}
				
				@Override
				public int getDomainDimension()
				{
					return face.getNormal().getDomainDimension();
				}
				
				@Override
				public CoordinateVector value(CoordinateVector pos)
				{
					return face.getNormal().value(pos).mul(-1);
				}
			};
		else return face.getNormal();
	}
	
	@Override
	public List<TPCell> refine(List<TPFace> refinedFaces)
	{
		throw new UnsupportedOperationException();
	}
	public String toString()
	{
		String ret = "";
		int subd = 0;
		for(int d = 0; d < getDimension(); d++)
		{
			ret = ret.concat("["+cell1Ds.get(subd).getStart()+", "+cell1Ds.get(subd++).getEnd()+ "]");
			if(d<getDimension()-1)
				ret = ret.concat("x");
		}
		return ret;
	}
	
	@Override
	public int compareTo(@NotNull TPCell o)
	{
		if(o.getDimension() < getDimension())
			return -1;
		if(o.getDimension() > getDimension())
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
		return ret;
	}
	
	@Override
	public boolean equals(Object obj)
	{
		if (obj instanceof TPCell)
			return 0 == compareTo((TPCell) obj);
		else
			return false;
	}
	
	@Override
	public void addEdge(TPEdge tpEdge)
	{
		
		if(edges.add(tpEdge))
			tpEdge.addCell(this);
	}
	@Override
	public TPCell getReferenceCell()
	{
		List<Cell1D> cells = new ArrayList<>(getDimension());
		for(int i = 0; i < getDimension(); i++)
			cells.add(new Cell1D(0,1));
		return new TPCell(cells);
	}
	
	
	@Override
	public CoordinateMatrix transformationGradientFromReferenceCell(CoordinateVector pos)
	{
		CoordinateMatrix ret = new CoordinateMatrix(getDimension(), getDimension());
		for(int i = 0; i < getDimension(); i++)
		{
			ret.set(cell1Ds.get(i).getEnd() - cell1Ds.get(i).getStart(), i, i);
		}
		return ret;
	}
	
	@Override
	public CoordinateMatrix transformationGradientToReferenceCell(CoordinateVector pos)
	{
		CoordinateMatrix ret = new CoordinateMatrix(getDimension(), getDimension());
		for(int i = 0; i < getDimension(); i++)
		{
			ret.set(1./(cell1Ds.get(i).getEnd() - cell1Ds.get(i).getStart()), i, i);
		}
		return ret;
	}
	
	@Override
	public CoordinateVector transformToReferenceCell(CoordinateVector pos)
	{
		CoordinateVector ret = new CoordinateVector(getDimension());
		for(int i = 0; i < getDimension(); i++)
		{
			ret.set((pos.at(i) - cell1Ds.get(i).getStart())/(cell1Ds.get(i).getEnd() - cell1Ds.get(i).getStart()), i);
		}
		return ret;
	}
	
	@Override
	public CoordinateVector transformFromReferenceCell(CoordinateVector pos)
	{
		CoordinateVector ret = new CoordinateVector(getDimension());
		for(int i = 0; i < getDimension(); i++)
		{
			ret.set(pos.at(i) *(cell1Ds.get(i).getEnd() - cell1Ds.get(i).getStart()) + cell1Ds.get(i).getStart(), i);
		}
		return ret;
	}
}



