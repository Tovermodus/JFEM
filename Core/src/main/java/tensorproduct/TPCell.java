package tensorproduct;

import basic.Cell;
import basic.Face;
import basic.ScalarShapeFunction;
import linalg.DoubleTensor;

import java.util.ArrayList;

public class TPCell extends Cell
{
	Cell1D cellx;
	Cell1D celly;
	private int polynomialDegree;
	public TPCell(double xStart, double yStart, double xEnd, double yEnd, int polynomialDegree)
	{
		super();
		cellx = new Cell1D(xStart, xEnd);
		celly = new Cell1D(yStart, yEnd);
		this.polynomialDegree = polynomialDegree;
	}
	public TPCell(Cell1D cellx,Cell1D celly, int polynomialDegree)
	{
		this.cellx = cellx;
		this.celly = celly;
		this.polynomialDegree = polynomialDegree;
	}

	@Override
	public void addFace(Face face)
	{
		this.faces.add(face);
		face.getCells().add(this);
	}

	@Override
	public void distributeFunctions(ArrayList<ScalarShapeFunction> globalshapeFunctions)
	{
		for(int i = 0; i <= polynomialDegree; i++)
		{
			for(int j = 0; j <= polynomialDegree; j++)
				shapeFunctions.add(new TPShapeFunction(new LagrangeBasisFunction1D(polynomialDegree,
					i, cellx), new LagrangeBasisFunction1D(polynomialDegree,j,celly),this));
		}
		for(ScalarShapeFunction shapeFunction:shapeFunctions)
		{
			for(Face face:faces)
			{
				face.getShapeFunctions().add(shapeFunction);
			}
			globalshapeFunctions.add(shapeFunction);
			shapeFunction.setGlobalIndex(globalshapeFunctions.size()-1);

		}
	}


	@Override
	public boolean isInCell(DoubleTensor pos)
	{
		return pos.x()>=cellx.getStart() && pos.x()<= cellx.getEnd() && pos.y()>=celly.getStart() && pos.y()<=celly.getEnd();
	}

	@Override
	public DoubleTensor center()
	{
		return DoubleTensor.vectorFromValues(0.5*(cellx.getStart() +cellx.getEnd()),0.5*(celly.getStart() +celly.getEnd()));
	}

	@Override
	public ArrayList<Cell> refine(ArrayList<Face> refinedFaces)
	{
		ArrayList<Cell> refinedCells = new ArrayList<>();
		TPCell cell1 = new TPCell(cellx.getStart(),celly.getStart(),cellx.center(),celly.center(),polynomialDegree);
		TPFace face1 = new TPFace(new Cell1D(cellx.getStart(),cellx.center()),celly.center(),1);
		TPFace face2 = new TPFace(new Cell1D(celly.getStart(),celly.center()),cellx.center(),0);
		TPCell cell2 = new TPCell(cellx.center(),celly.getStart(),cellx.getEnd(),celly.center(),polynomialDegree);
		TPCell cell3 = new TPCell(cellx.getStart(),celly.center(),cellx.center(),celly.getEnd(),polynomialDegree);
		TPCell cell4 = new TPCell(cellx.center(),celly.center(),cellx.getEnd(),celly.getEnd(),polynomialDegree);
		TPFace face3 = new TPFace(new Cell1D(cellx.center(),cellx.getEnd()),celly.center(),1);
		TPFace face4 = new TPFace(new Cell1D(celly.center(),celly.getEnd()),cellx.center(),0);
		cell1.faces.add(face1);
		cell1.faces.add(face2);
		cell2.faces.add(face2);
		cell2.faces.add(face3);
		cell3.faces.add(face1);
		cell3.faces.add(face4);
		cell4.faces.add(face3);
		cell4.faces.add(face4);
		System.out.println("MAKE NICER TPCELL");
		face1.getCells().add(cell1);
		face1.getCells().add(cell3);
		face2.getCells().add(cell1);
		face2.getCells().add(cell2);
		face3.getCells().add(cell2);
		face3.getCells().add(cell4);
		face4.getCells().add(cell3);
		face4.getCells().add(cell4);
		face1.setBoundaryFace(false);
		face2.setBoundaryFace(false);
		face3.setBoundaryFace(false);
		face4.setBoundaryFace(false);
		refinedCells.add(cell1);
		refinedCells.add(cell2);
		refinedCells.add(cell3);
		refinedCells.add(cell4);
		refinedFaces.add(face1);
		refinedFaces.add(face2);
		refinedFaces.add(face3);
		refinedFaces.add(face4);
		setRefined(true);
		return refinedCells;
	}
	public void print()
	{
		System.out.println("["+cellx.getStart()+","+cellx.getEnd()+"]×["+celly.getStart()+","+celly.getEnd()+"]");
	}


}



