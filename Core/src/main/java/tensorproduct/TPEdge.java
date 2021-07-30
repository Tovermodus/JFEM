package tensorproduct;

import basic.DoubleCompare;
import basic.EdgeWithReferenceEdge;
import basic.PerformanceArguments;
import basic.VectorFunction;
import linalg.CoordinateComparator;
import linalg.CoordinateVector;
import org.jetbrains.annotations.NotNull;

import java.util.Set;
import java.util.TreeSet;

public class TPEdge implements EdgeWithReferenceEdge<TPCell, TPFace, TPEdge>
{
	double[] otherCoordinates;
	Cell1D cell;
	int tangentialDimension;
	private final Set<TPCell> cells;
	private final Set<TPFace> faces;
	private final VectorFunction tangent;
	boolean isBoundaryFace = false;
	
	public TPEdge(Cell1D cell, double [] otherCoordinates, int tangentialDimension)
	{
		if(otherCoordinates.length != 2)
			throw new IllegalArgumentException("onlly in 3D");
		this.tangentialDimension = tangentialDimension;
		this.cell = new Cell1D(cell);
		this.otherCoordinates = otherCoordinates.clone();
		this.cells = new TreeSet<>();
		this.faces = new TreeSet<>();
		this.tangent = new VectorFunction()
		{
			@Override
			public int getRangeDimension()
			{
				return 3;
			}
			
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
	public static TPEdge createEdgeFromFace(TPFace face, int eliminatedDirection,  boolean leftSide)
	{
		int retainedDirection = 0;
		for(int d = 0; d < 3; d ++)
			if(d != face.getFlatDimension() && d != eliminatedDirection)
				retainedDirection = d;
		
		Cell1D retainedCell = face.getCell1Ds().get(retainedDirection>eliminatedDirection?1:0);
		Cell1D eliminatedCell = face.getCell1Ds().get(retainedDirection>eliminatedDirection?0:1);
		
		double[] otherCoordinates = new double[]{-1,-1};
		if(eliminatedDirection < face.getFlatDimension())
		{
			if(leftSide)
				otherCoordinates[0] = eliminatedCell.getStart();
			else
				otherCoordinates[0] = eliminatedCell.getEnd();
			otherCoordinates[1] = face.getOtherCoordinate();
		}
		if(eliminatedDirection > face.getFlatDimension())
		{
			if(leftSide)
				otherCoordinates[1] = eliminatedCell.getStart();
			else
				otherCoordinates[1] = eliminatedCell.getEnd();
			otherCoordinates[0] = face.getOtherCoordinate();
		}
		TPEdge e = new TPEdge(retainedCell, otherCoordinates, retainedDirection);
		e.addFace(face);
		return e;
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
		if(face.isBoundaryFace())
			this.setBoundaryFace(true);
		if(faces.add(face))
			face.addEdge(this);
		for(TPCell c: face.getCells())
			for(TPFace f: c.getFaces())
				if(f.isOnFace(center()))
					if(!faces.contains(f))
					{
						addFace(f);
					}
		for(TPCell c: face.getCells())
			addCell(c);
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
				if (!DoubleCompare.almostEqual(otherCoordinates[subd++],  pos.at(d)))
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
		if(o.otherCoordinates[0] < otherCoordinates[0]-PerformanceArguments.getInstance().doubleTolerance)
			return -1;
		if(o.otherCoordinates[0] > otherCoordinates[0]+PerformanceArguments.getInstance().doubleTolerance)
			return 1;
		if(o.otherCoordinates[1] < otherCoordinates[1]-PerformanceArguments.getInstance().doubleTolerance)
			return -1;
		if(o.otherCoordinates[1] > otherCoordinates[1]+PerformanceArguments.getInstance().doubleTolerance)
			return 1;
		return CoordinateComparator.comp(center().getEntries(), o.center().getEntries());
	}
	
	@Override
	public int hashCode()
	{
		int ret = 0;
		ret += Math.pow(7, 1) * cell.center();
		ret += Math.pow(7, 2) * cell.getStart();
		ret += Math.pow(7, 3) * cell.getEnd();
		ret -= 141*otherCoordinates[0];
		ret -= 141*141*otherCoordinates[1];
		return ret;
	}
	
	@Override
	public boolean equals(Object obj)
	{
		if (obj instanceof TPEdge)
			return 0 == compareTo((TPEdge) obj);
		else
			return false;
	}
	@Override
	public void setBoundaryFace(boolean boundaryFace)
	{
		isBoundaryFace = boundaryFace;
	}
	
	@Override
	public boolean isBoundaryEdge()
	{
		return isBoundaryFace;
	}
	
	public Cell1D getCell()
	{
		return cell;
	}
	
	public int getTangentialDimension()
	{
		return tangentialDimension;
	}
	
	@Override
	public TPEdge getReferenceEdge()
	{
		double [] othercoordinates = {0,0};
		return new TPEdge(new Cell1D(0,1), othercoordinates, tangentialDimension);
	}
}
