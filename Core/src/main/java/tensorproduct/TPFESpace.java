package tensorproduct;

import basic.Assembleable;
import basic.Cell;
import basic.FESpace;

import java.util.ArrayList;

public class TPFESpace extends FESpace implements Assembleable
{
	/*
	TODO: Provide interface for fast TPIntegrals
	 */
	ArrayList<Cell1D> cellsX;
	ArrayList<Cell1D> cellsY;

	Cell[][] cellsArray;
	double xStart;
	double yStart;
	double xEnd;
	double yEnd;
	int numberXCells;
	int numberYCells;
	public TPFESpace(double xStart, double yStart, double xEnd, double yEnd, int numberXCells, int numberYCells,
	                 int polynomialDegree)
	{
		super();
		this.xStart = xStart;
		this.xEnd = xEnd;
		this.yStart = yStart;
		this.yEnd = yEnd;
		this.numberXCells = numberXCells;
		this.numberYCells = numberYCells;
		cellsX = new ArrayList<>();
		cellsY = new ArrayList<>();
		assembleCells(polynomialDegree);
		assembleFunctions();
	}
	@Override
	public void assembleFunctions()
	{
		for(Cell cell: cells)
		{
			cell.distributeFunctions(shapeFunctions);
		}
	}
	@Override
	public void assembleCells(int polynomialDegree)
	{
		cellsArray = new Cell[numberXCells][numberYCells];
		double hx =  (xEnd - xStart)/numberXCells;
		double hy =  (yEnd - yStart)/numberYCells;
		for(int i = 0; i < numberXCells; i++)
		{
			cellsX.add(new Cell1D(xStart + i * hx, xStart + (i + 1) * hx));
		}
		for(int i = 0; i < numberYCells; i++)
			cellsY.add(new Cell1D(yStart+i*hy,yStart+(i+1)*hy));
		for(int i = 0; i < numberXCells; i++)
		{
			System.out.println("assemble cells: "+(int)(1.0*i/(numberXCells)*100)+"%");
			for(int j = 0; j < numberYCells; j++)
			{
				cellsArray[i][j] = new TPCell(cellsX.get(i), cellsY.get(j),polynomialDegree);
				cells.add(cellsArray[i][j]);
			}
		}
		for(int i = 0; i < numberXCells+1; i++)
		{
			System.out.println("assemble faces: "+(int)((1.0*i+1)/(numberXCells+1)*100)+"%");
			for(int j = 0; j < numberYCells+1; j++)
			{
				if(j < numberYCells)
				{
					TPFace xFace = new TPFace(cellsY.get(j),xStart+i*hx,0); //normal direction is
					// x axis
					faces.add(xFace);
					if(i == 0)
					{
						xFace.setBoundaryFace(true);
						cellsArray[i][j].addFace(xFace);
					}
					else if(i == numberXCells)
					{
						xFace.setBoundaryFace(true);
						cellsArray[i-1][j].addFace(xFace);
					}
					else
					{
						cellsArray[i][j].addFace(xFace);
						cellsArray[i-1][j].addFace(xFace);
					}
				}
				if(i < numberXCells)
				{
					TPFace yFace = new TPFace(cellsX.get(i), yStart + j * hy, 1);
					faces.add(yFace);
					if (j == 0)
					{
						yFace.setBoundaryFace(true);
						cellsArray[i][j].addFace(yFace);
					} else if (j == numberYCells)
					{
						yFace.setBoundaryFace(true);
						cellsArray[i][j - 1].addFace(yFace);
					} else
					{
						cellsArray[i][j].addFace(yFace);
						cellsArray[i][j - 1].addFace(yFace);
					}
				}
			}
		}
	}
}
