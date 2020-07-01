package tensorproduct;

import basic.Face;
import basic.TensorFunction;
import com.google.common.collect.Multimap;
import linalg.DoubleTensor;
import basic.Cell;
import java.util.ArrayList;

public class TPFace extends Face
{
	private Cell1D cell1d;
	private double otherCoordinate;
	private int normaldirection;
	public TPFace(Cell1D cell1d, double otherCoordinate, int normaldirection)
	{
		super(new TensorFunction()
		{
			@Override
			public DoubleTensor value(DoubleTensor pos)
			{
				if(normaldirection == 0)
					return DoubleTensor.vectorFromValues(1,0);
				else
					return DoubleTensor.vectorFromValues(0,1);
			}

			@Override
			public DoubleTensor derivative(DoubleTensor pos)
			{
				return DoubleTensor.squareMatrixFromValues(0,0,0,0);
			}
		});
		this.cell1d = cell1d;
		this.otherCoordinate = otherCoordinate;
		this.normaldirection = normaldirection;
	}
	
	public Cell1D getCell1d()
	{
		return cell1d;
	}
	
	public int getNormaldirection()
	{
		return normaldirection;
	}
	
	public double getOtherCoordinate()
	{
		return otherCoordinate;
	}
	
	public Cell getUpStreamCell(DoubleTensor pos, DoubleTensor direction)
	{
		if(normal.value(pos).inner(direction)>0)
			return getNormalUpstreamCell(pos);
		else
			return getNormalDownstreamCell(pos);
	}
	public Cell getDownStreamCell(DoubleTensor pos, DoubleTensor direction)
	{
		if(normal.value(pos).inner(direction)<0)
			return getNormalUpstreamCell(pos);
		else
			return getNormalDownstreamCell(pos);
	}
	@Override
	public Cell getNormalUpstreamCell(DoubleTensor pos)
	{
		for(Cell cell:cells)
		{
			if(cell.isInCell(pos))
			{
				if (cell.center().at(normaldirection) < otherCoordinate)
				{
					return cell;
				}
			}
		}
		return null;
	}

	@Override
	public Cell getNormalDownstreamCell(DoubleTensor pos)
	{
		for(Cell cell:cells)
		{

			if(cell.isInCell(pos))
			{
				if (cell.center().at(normaldirection) > otherCoordinate)
					return cell;
			}
		}
		return null;
	}

	@Override
	public DoubleTensor center()
	{
		if(normaldirection == 0)
			return DoubleTensor.vectorFromValues(otherCoordinate,cell1d.center());
		else
			return DoubleTensor.vectorFromValues(cell1d.center(),otherCoordinate);

	}

	@Override
	public boolean isOnFace(DoubleTensor pos)
	{
		return pos.at(normaldirection) == otherCoordinate && cell1d.isInCell(pos.at(1-normaldirection));
	}

	@Override
	public ArrayList<Face> refine(Multimap<Cell, Cell> cellMap)
	{
		assert(cells.size()<2);
		ArrayList<Face> refinedFaces = new ArrayList<>();
		TPFace face1 = new TPFace(new Cell1D(cell1d.getStart(),cell1d.center()),otherCoordinate,normaldirection);
		TPFace face2 = new TPFace(new Cell1D(cell1d.center(), cell1d.getEnd()),otherCoordinate,normaldirection);
		for(Cell cell:cells)
		{
			for (Cell refinedCell : cellMap.get(cell))
			{
				if(refinedCell.isInCell(face1.center()))
				{
					refinedCell.getFaces().add(face1);
					face1.cells.add(refinedCell);
				}
				if(refinedCell.isInCell(face2.center()))
				{
					refinedCell.getFaces().add(face2);
					face2.cells.add(refinedCell);
				}
			}
		}
		face1.setBoundaryFace(isBoundaryFace());
		face2.setBoundaryFace(isBoundaryFace());
		refinedFaces.add(face1);
		refinedFaces.add(face2);
		return refinedFaces;

	}
    public void print()
	{
		if(this.normaldirection == 0)
		{
			System.out.println(otherCoordinate+"×["+ cell1d.getStart() +","+ cell1d.getEnd() +"]   " +cells.size());
		}
		else
		{
			System.out.println("["+ cell1d.getStart() +","+ cell1d.getEnd() +"]×"+otherCoordinate+"   "+cells.size());
		}
	       
	}



}
