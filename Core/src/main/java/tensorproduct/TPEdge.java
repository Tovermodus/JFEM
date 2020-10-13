package tensorproduct;

import basic.Edge;
import basic.VectorFunction;
import linalg.CoordinateComparator;
import linalg.CoordinateVector;
import org.jetbrains.annotations.NotNull;

import java.util.Set;
import java.util.TreeSet;

public class TPEdge implements Edge<TPCell, TPFace, TPEdge>
{
	double[] otherCoordinates;
	Cell1D cell;
	int tangentialDimension;
	private Set<TPCell> cells;
	private Set<TPFace> faces;
	private final VectorFunction tangent;
	
	public TPEdge(Cell1D cell, double [] otherCoordinates, int tangentialDimension)
	{
		if(otherCoordinates.length != 2)
			throw new IllegalArgumentException("onlly in 3D");
		this.tangentialDimension = tangentialDimension;
		this.cell = cell;
		this.otherCoordinates = otherCoordinates;
		this.cells = new TreeSet<>();
		this.faces = new TreeSet<>();
		this.tangent = new VectorFunction()
		{
			@Override
			public int getDomainDimension()
			{
				return 3;
			}
			
			@Override
			public CoordinateVector value(CoordinateVector pos)
			{
				return CoordinateVector.getUnitVector(3,tangentialDimension);
			}
		};
	}
	
	@Override
	public void addCell(TPCell cell)
	{
		if(cells.add(cell))
			cell.addEdge(this);
	}
	@Override
	public void addFace(TPFace face)
	{
		if(faces.add(face))
			face.addEdge(this);
	}
	@Override
	public int getDimension()
	{
		return 3;
	}
	
	@Override
	public Set<TPCell> getCells()
	{
		return cells;
	}
	
	@Override
	public Set<TPFace> getFaces()
	{
		return faces;
	}
	
	@Override
	public VectorFunction getTangent()
	{
		return tangent;
	}
	
	@Override
	public CoordinateVector center()
	{
		
		CoordinateVector ret = new CoordinateVector(3);
		int subd = 0;
		for(int d = 0; d < getDimension(); d++)
		{
			if(d == tangentialDimension)
				ret.set(cell.center(), d);
			else
				ret.set(otherCoordinates[subd++],d);
		}
		return ret;
	}
	
	@Override
	public boolean isOnEdge(CoordinateVector pos)
	{
		int subd = 0;
		for(int d = 0; d < getDimension()+1; d++)
		{
			if(d == tangentialDimension)
			{
				if (!cell.isInCell(pos.at(tangentialDimension)))
					return false;
			}
			else
			{
				if (otherCoordinates[subd++] != pos.at(d))
					return false;
			}
		}
		return true;
		
	}
	
	@Override
	public int compareTo(@NotNull TPEdge o)
	{
		if(o.tangentialDimension < tangentialDimension)
			return -1;
		if(o.tangentialDimension > tangentialDimension)
			return 1;
		if(o.otherCoordinates[0] < otherCoordinates[0])
			return -1;
		if(o.otherCoordinates[0] > otherCoordinates[0])
			return 1;
		if(o.otherCoordinates[1] < otherCoordinates[1])
			return -1;
		if(o.otherCoordinates[1] > otherCoordinates[1])
			return 1;
		return CoordinateComparator.comp(center().getEntries(), o.center().getEntries());
	}
	
	public Cell1D getCell()
	{
		return cell;
	}
	
	public int getTangentialDimension()
	{
		return tangentialDimension;
	}
}
