package tensorproduct.geometry;

import basic.*;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableSet;
import linalg.CoordinateComparator;
import linalg.CoordinateMatrix;
import linalg.CoordinateVector;
import org.jetbrains.annotations.NotNull;

import java.util.*;

public class TPCell implements CellWithReferenceCell<TPCell, TPFace>
{
	final ImmutableList<Cell1D> cell1Ds;
	Set<TPFace> faces;
	boolean refined;
	
	TPCell(List<Cell1D> cell1Ds)
	{
		this.cell1Ds = ImmutableList.copyOf(cell1Ds);
		this.faces = new HashSet<>(2*getDimension());
		this.refined = false;
	}
	
	@Override
	public int getDimension()
	{
		return cell1Ds.size();
	}
	
	public Cell1D getComponentCell(int dim)
	{
		return cell1Ds.get(dim);
	}
	public ImmutableList<Cell1D> getComponentCells()
	{
		return cell1Ds;
	}
	
	@Override
	public ImmutableSet<TPFace> getFaces()
	{
		return ImmutableSet.copyOf(faces);
	}
	
	@Override
	public boolean isInCell(CoordinateVector pos)
	{
		for (int d = 0; d < cell1Ds.size(); d++)
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
		for (int i = 0; i < ret.getLength(); i++)
			ret.set(cell1Ds.get(i).center(), i);
		return ret;
	}
	
	@Override
	public VectorFunction getOuterNormal(TPFace face)
	{
		if (PerformanceArguments.getInstance().executeChecks)
			if (!faces.contains(face))
				throw new IllegalArgumentException("face does not belong to cell");
		boolean invertNormal = (center().sub(face.center())).inner(face.getNormal().value(face.center())) > 0;
		if (invertNormal)
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
		for (int d = 0; d < getDimension(); d++)
		{
			ret = ret.concat("[" + cell1Ds.get(subd).getStart() + ", " + cell1Ds.get(subd++).getEnd() + "]");
			if (d < getDimension() - 1)
				ret = ret.concat("x");
		}
		return ret;
	}
	
	@Override
	public int compareTo(@NotNull TPCell o)
	{
		if (o.getDimension() < getDimension())
			return -1;
		if (o.getDimension() > getDimension())
			return 1;
		for (int i = 0; i < cell1Ds.size(); i++)
		{
			if(!DoubleCompare.almostEqual(getComponentCell(i).getStart(),
				o.getComponentCell(i).getStart()))
				return Double.compare(getComponentCell(i).getStart(),
					o.getComponentCell(i).getStart());
			if(!DoubleCompare.almostEqual(getComponentCell(i).getEnd(),
				o.getComponentCell(i).getEnd()))
				return Double.compare(getComponentCell(i).getEnd(),
					o.getComponentCell(i).getEnd());
		}
		return CoordinateComparator.comp(center().getEntries(), o.center().getEntries());
	}
	
	@Override
	public int hashCode()
	{
		int ret = 0;
		for (int i = 0; i < cell1Ds.size(); i++)
		{
			ret += 175628373 *3*i* DoubleCompare.doubleHash(cell1Ds.get(i).center());
			ret += 175628373*(3*i+1)* DoubleCompare.doubleHash(cell1Ds.get(i).getStart());
			ret += 175628373 *(3*i+2)* DoubleCompare.doubleHash(cell1Ds.get(i).getEnd());
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
	public TPCell getReferenceCell()
	{
		List<Cell1D> cells = new ArrayList<>(getDimension());
		for (int i = 0; i < getDimension(); i++)
			cells.add(new Cell1D(0, 1));
		return new TPCell(cells);
	}
	
	
	@Override
	public CoordinateMatrix transformationGradientFromReferenceCell(CoordinateVector pos)
	{
		CoordinateMatrix ret = new CoordinateMatrix(getDimension(), getDimension());
		for (int i = 0; i < getDimension(); i++)
		{
			ret.set(cell1Ds.get(i).getEnd() - cell1Ds.get(i).getStart(), i, i);
		}
		return ret;
	}
	
	@Override
	public CoordinateMatrix transformationGradientToReferenceCell(CoordinateVector pos)
	{
		CoordinateMatrix ret = new CoordinateMatrix(getDimension(), getDimension());
		for (int i = 0; i < getDimension(); i++)
		{
			ret.set(1. / (cell1Ds.get(i).getEnd() - cell1Ds.get(i).getStart()), i, i);
		}
		return ret;
	}
	
	@Override
	public CoordinateVector transformToReferenceCell(CoordinateVector pos)
	{
		CoordinateVector ret = new CoordinateVector(getDimension());
		for (int i = 0; i < getDimension(); i++)
		{
			ret.set((pos.at(i) - cell1Ds.get(i).getStart()) / (cell1Ds.get(i).getEnd() - cell1Ds.get(i).getStart()), i);
		}
		return ret;
	}
	
	@Override
	public CoordinateVector transformFromReferenceCell(CoordinateVector pos)
	{
		CoordinateVector ret = new CoordinateVector(getDimension());
		for (int i = 0; i < getDimension(); i++)
		{
			ret.set(pos.at(i) * (cell1Ds.get(i).getEnd() - cell1Ds.get(i).getStart()) + cell1Ds.get(i).getStart(), i);
		}
		return ret;
	}
}



