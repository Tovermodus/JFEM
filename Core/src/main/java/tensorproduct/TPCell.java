package tensorproduct;

import basic.Cell;
import basic.ShapeFunction;
import linalg.CoordinateComparator;
import linalg.CoordinateVector;
import org.jetbrains.annotations.NotNull;

import java.util.List;
import java.util.Set;
import java.util.TreeSet;

public class TPCell implements Cell<TPCell, TPFace>
{
	List<Cell1D> cell1Ds;
	Set<TPFace> faces;
	boolean refined;
	
	public TPCell(List<Cell1D> cell1Ds)
	{
		this.cell1Ds = cell1Ds;
		this.faces = new TreeSet<>();
		this.refined = false;
	}
	@Override
	public int getDimension()
	{
		return cell1Ds.size();
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
	
	//	@Override
//	public  List<TPCell> refine(List<TPFace> refinedFaces)
//	{
//		List<Cell<TPFace,TPShapeFunction>> refinedCells = new ArrayList<>();
//		TPCell cell1 = new TPCell(cellx.getStart(),celly.getStart(),cellx.center(),celly.center(),
//			polynomialDegree);
//		TPFace face1 = new TPFace(new Cell1D(cellx.getStart(),cellx.center()),celly.center(),1);
//		TPFace face2 = new TPFace(new Cell1D(celly.getStart(),celly.center()),cellx.center(),0);
//		TPCell cell2 = new TPCell(cellx.center(),celly.getStart(),cellx.getEnd(),celly.center(),
//			polynomialDegree);
//		TPCell cell3 = new TPCell(cellx.getStart(),celly.center(),cellx.center(),celly.getEnd(),
//			polynomialDegree);
//		TPCell cell4 = new TPCell(cellx.center(),celly.center(),cellx.getEnd(),celly.getEnd(),
//			polynomialDegree);
//		TPFace face3 = new TPFace(new Cell1D(cellx.center(),cellx.getEnd()),celly.center(),1);
//		TPFace face4 = new TPFace(new Cell1D(celly.center(),celly.getEnd()),cellx.center(),0);
//		cell1.faces.add(face1);
//		cell1.faces.add(face2);
//		cell2.faces.add(face2);
//		cell2.faces.add(face3);
//		cell3.faces.add(face1);
//		cell3.faces.add(face4);
//		cell4.faces.add(face3);
//		cell4.faces.add(face4);
//		System.out.println("MAKE NICER TPCELL");
//		face1.getCells().add(cell1);
//		face1.getCells().add(cell3);
//		face2.getCells().add(cell1);
//		face2.getCells().add(cell2);
//		face3.getCells().add(cell2);
//		face3.getCells().add(cell4);
//		face4.getCells().add(cell3);
//		face4.getCells().add(cell4);
//		face1.setBoundaryFace(false);
//		face2.setBoundaryFace(false);
//		face3.setBoundaryFace(false);
//		face4.setBoundaryFace(false);
//		refinedCells.add(cell1);
//		refinedCells.add(cell2);
//		refinedCells.add(cell3);
//		refinedCells.add(cell4);
//		refinedFaces.add(face1);
//		refinedFaces.add(face2);
//		refinedFaces.add(face3);
//		refinedFaces.add(face4);
//		setRefined(true);
//		return refinedCells;
//	}
//


}


