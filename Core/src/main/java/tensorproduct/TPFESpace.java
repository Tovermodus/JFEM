package tensorproduct;

import basic.*;
import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;
import com.google.common.collect.Multimap;
import com.google.common.primitives.Doubles;
import com.google.common.primitives.Ints;
import linalg.*;
import linalg.Vector;

import java.awt.image.AreaAveragingScaleFilter;
import java.util.*;
import java.util.stream.Collectors;

/*
Workflow of this finite element is: first generate 1D cells, then
 */
public class TPFESpace implements MatrixFESpace<TPCell<TPShapeFunction>,
	TPFace<TPShapeFunction>,
	TPShapeFunction,Double, Vector,Matrix,TPFESpace>,
	Assembleable
{
	List<List<Double>> coordinates1D;
	List<List<LagrangeBasisFunction1D>> basisFunctions;
	List<List<Cell1D>> cells1D;
	List<TPCell<TPShapeFunction>> cells;
	List<TPFace<TPShapeFunction>> faces;
	Map<List<Integer>, TPCell<TPShapeFunction>> lexicographicCellNumbers;
	List<TPShapeFunction> shapeFunctions;
	DenseMatrix systemMatrix;
	DenseVector rhs;
	public TPFESpace(CoordinateVector startCoordinates, CoordinateVector endCoordinates,
	                 List<Integer> cellsPerDimension, int polynomialDegree)
	{
		if(startCoordinates.getLength() != endCoordinates.getLength()|| startCoordinates.getLength() != cellsPerDimension.size())
			throw new IllegalArgumentException();
		cells1D = new ArrayList<>();
		coordinates1D = new ArrayList<>();
		QuadratureRule1D quad;
		if(polynomialDegree < 3)
			quad = QuadratureRule1D.Gauss3;
		else
			quad = QuadratureRule1D.Gauss5;
		for (int i = 0; i < startCoordinates.getLength(); i++)
		{
			List<Cell1D> cellsForDirection = new ArrayList<>();
			List<Double> coordinatesForDirection = new ArrayList<>();
			double le = (endCoordinates.at(i) - startCoordinates.at(i))/cellsPerDimension.get(i);
			for (int j = 0; j < cellsPerDimension.get(i); j++)
			{
				coordinatesForDirection.add(startCoordinates.at(i) + le*j);
				cellsForDirection.add(new Cell1D(startCoordinates.at(i) + le*j,
					startCoordinates.at(i) + le*(j+1),quad, cellsForDirection.size()));
			}
			coordinatesForDirection.add(endCoordinates.at(i));
			cells1D.add(cellsForDirection);
			coordinates1D.add(coordinatesForDirection);
		}
		
	}
	@Override
	public void assembleFunctions(int polynomialDegree)
	{
		shapeFunctions = new ArrayList<>();
		basisFunctions = new ArrayList<>();
		for (int i = 0; i < coordinates1D.size(); i++)
		{
			List<LagrangeBasisFunction1D> functionsForDirection = new ArrayList<>();
			for (int j = 0; j < cells1D.get(i).size(); j++)
			{
				for (int k = 0; k < polynomialDegree + 1; k++)
				{
					functionsForDirection.add(new LagrangeBasisFunction1D(polynomialDegree, k,
						cells1D.get(i).get(j)));
				}
			}
			basisFunctions.add(functionsForDirection);
		}
		for(List<LagrangeBasisFunction1D> functionCombination: Lists.cartesianProduct(basisFunctions))
		{
			int [] lexicographicCellIndex =
				functionCombination.stream().mapToInt(LagrangeBasisFunction1D::getCellIndexInDimension).toArray();
			System.out.println(Arrays.deepToString(Ints.asList(lexicographicCellIndex).toArray()));
			shapeFunctions.add(new TPShapeFunction(lexicographicCellNumbers.get(Ints.asList(lexicographicCellIndex)),
				functionCombination));
			Iterables.getLast(shapeFunctions).setGlobalIndex(shapeFunctions.size()-1);
			System.out.println(Arrays.deepToString(shapeFunctions.toArray()));
		}
	}
	@Override
	public void assembleCells()
	{
		cells = new ArrayList<>();
		faces = new ArrayList<>();
		lexicographicCellNumbers = new HashMap<>();
		for(List<Cell1D> cellCombination: Lists.cartesianProduct(cells1D))
		{
			cells.add(new TPCell<>(cellCombination));
			lexicographicCellNumbers.put(cellCombination.stream().map(Cell1D::getIndexInDimension).collect(Collectors.toList()),
				Iterables.getLast(cells));
		}
		for(int i = 0; i < cells1D.size(); i++)
		{
			List<List<Cell1D>> cells1Dcopy = new ArrayList<>(cells1D);
			cells1Dcopy.remove(i);
			for(int j = 0; j < coordinates1D.get(i).size(); j++)
			{
				for(List<Cell1D> cellCombination: Lists.cartesianProduct(cells1Dcopy))
				{
					faces.add(new TPFace<>(cellCombination,i,coordinates1D.get(i).get(j),
						(j == 0 || j == coordinates1D.get(i).size()-1)));
					int [] lexicographicIndex = new int[coordinates1D.size()];
					int[] subLexIndex =
						cellCombination.stream().mapToInt(Cell1D::getIndexInDimension).toArray();
					int subd = 0;
					for(int k = 0; k < coordinates1D.size(); k++)
					{
						if(k == i)
							lexicographicIndex[k] = j-1;
						else
							lexicographicIndex[k] = subLexIndex[subd++];
					}
					if(j!=0)
					{
						lexicographicCellNumbers.get(Ints.asList(lexicographicIndex)).addFace(Iterables.getLast(faces));
					}
					if(j != cells1D.get(i).size())
					{
						lexicographicIndex[i] = j;
						lexicographicCellNumbers.get(Ints.asList(lexicographicIndex)).addFace(Iterables.getLast(faces));
					}
					
				}
			}
		}
		for(TPCell<TPShapeFunction> cell: cells)
			System.out.println(Arrays.deepToString(cell.getFaces().toArray()));
	}
	
	@Override
	public void initializeSystemMatrix()
	{
		systemMatrix = new DenseMatrix(shapeFunctions.size(),shapeFunctions.size());
	}
	
	@Override
	public void initializeRhs()
	{
		rhs = new DenseVector(shapeFunctions.size());
	}
	
	@Override
	public Vector getRhs()
	{
		return rhs;
	}
	
	@Override
	public Matrix getSystemMatrix()
	{
		return systemMatrix;
	}
	
	@Override
	public List<TPCell<TPShapeFunction>> getCells()
	{
		return cells;
	}
	
	@Override
	public List<TPShapeFunction> getShapeFunctions()
	{
		return shapeFunctions;
	}
	
	@Override
	public List<TPFace<TPShapeFunction>> getFaces()
	{
		return faces;
	}
	
	@Override
	public List<CoordinateVector> generatePlotPoints(double resolution)
	{
		int n = (int)(1./(1.-resolution));
		List<List<Double>> plotCoordinates1D = new ArrayList<>();
		for (int i = 0; i < coordinates1D.size(); i++)
		{
			List<Double> coordinatesForDirection = new ArrayList<>();
			double le =
				(Iterables.getLast(coordinates1D.get(i)) -coordinates1D.get(i).get(0))/n;
			for (int j = 0; j < n; j++)
			{
				coordinatesForDirection.add(coordinates1D.get(i).get(0) + le*j);
			}
			coordinatesForDirection.add(Iterables.getLast(coordinates1D.get(i)));
			plotCoordinates1D.add(coordinatesForDirection);
		}
		return Lists.cartesianProduct(plotCoordinates1D).stream().map(Doubles::toArray).map(CoordinateVector::fromValues).collect(Collectors.toList());
	}
	
	@Override
	public TPFESpace refine(Multimap<TPCell<TPShapeFunction>,TPCell<TPShapeFunction>> cellRefinedCellMapping,
	                        Multimap<TPFace<TPShapeFunction>,TPFace<TPShapeFunction>> faceRefinedFaceMapping)
	{
		return null;
	}


//	/*
//	TODO: Provide interface for fast TPIntegrals
//	 */
//
//	private ArrayList<TPCell> cells;
//	private ArrayList<TPFace> faces;
//	private ArrayList<TPShapeFunction> shapeFunctions;
//	private DoubleTensor systemMatrix;
//	private DoubleTensor rhs;
//	private final ArrayList<Cell1D> cellsX;
//	private final ArrayList<Cell1D> cellsY;
//	private final ArrayList<LagrangeBasisFunction1D> functionsX;
//	private final ArrayList<LagrangeBasisFunction1D> functionsY;
//
//	private TPCell[][] cellsArray;
//	private TPShapeFunction[][] functionsArray;
//	private final double xStart;
//	private final double yStart;
//	private final double xEnd;
//	private final double yEnd;
//	private final int numberXCells;
//	private final int numberYCells;
//	private int polynomialDegree;
//	public TPFESpace(double xStart, double yStart, double xEnd, double yEnd, int numberXCells, int numberYCells,
//	                 int polynomialDegree)
//	{
//		super();
//		this.xStart = xStart;
//		this.xEnd = xEnd;
//		this.yStart = yStart;
//		this.yEnd = yEnd;
//		this.numberXCells = numberXCells;
//		this.numberYCells = numberYCells;
//		cellsX = new ArrayList<>();
//		cellsY = new ArrayList<>();
//		functionsX = new ArrayList<>();
//		functionsY = new ArrayList<>();
//		assembleCells(polynomialDegree);
//		assembleFunctions();
//	}
//	@Override
//	public void assembleFunctions()
//	{
//		shapeFunctions = new ArrayList<>();
//		for(Cell1D cellX:getCellsX())
//			for(int j = 0; j < polynomialDegree; j++)
//				functionsX.add(new LagrangeBasisFunction1D(polynomialDegree, j, cellX));
//		for(Cell1D cellY:getCellsY())
//			for(int j = 0; j < polynomialDegree; j++)
//				functionsY.add(new LagrangeBasisFunction1D(polynomialDegree, j, cellY));
//
//		functionsArray = new TPShapeFunction[functionsX.size()][functionsY.size()];
//		for(int i = 0; i < functionsX.size(); i++)
//			for(int j = 0; j < functionsY.size(); j++)
//			{
//				TPCell cell = cellsArray[i/polynomialDegree][j/polynomialDegree];
//				functionsArray[i][j] = new TPShapeFunction(functionsX.get(i), functionsY.get(j),cell);
//				functionsArray[i][j].setGlobalIndex(shapeFunctions.size());
//				shapeFunctions.add(functionsArray[i][j]);
//			}
//	}
//	@Override
//	public void assembleCells(int polynomialDegree)
//	{
//		cellsArray = new TPCell[getNumberXCells()][getNumberYCells()];
//		double hx =  (getxEnd() - getxStart())/ getNumberXCells();
//		double hy =  (getyEnd() - getyStart())/ getNumberYCells();
//		for(int i = 0; i < getNumberXCells(); i++)
//		{
//			getCellsX().add(new Cell1D(getxStart() + i * hx, getxStart() + (i + 1) * hx));
//		}
//		for(int i = 0; i < getNumberYCells(); i++)
//			getCellsY().add(new Cell1D(getyStart() +i*hy, getyStart() +(i+1)*hy));
//		for(int i = 0; i < getNumberXCells(); i++)
//		{
//			System.out.println("assemble cells: "+(int)(1.0*i/(getNumberXCells())*100)+"%");
//			for(int j = 0; j < getNumberYCells(); j++)
//			{
//				getCellsArray()[i][j] = new TPCell(getCellsX().get(i), getCellsY().get(j),polynomialDegree);
//				getCells().add(getCellsArray()[i][j]);
//			}
//		}
//		for(int i = 0; i < getNumberXCells() +1; i++)
//		{
//			System.out.println("assemble faces: "+(int)((1.0*i+1)/(getNumberXCells() +1)*100)+"%");
//			for(int j = 0; j < getNumberYCells() +1; j++)
//			{
//				if(j < getNumberYCells())
//				{
//					TPFace xFace = new TPFace(getCellsY().get(j), getxStart() +i*hx,0); //normal direction is
//					// x axis
//					getFaces().add(xFace);
//					if(i == 0)
//					{
//						xFace.setBoundaryFace(true);
//						getCellsArray()[i][j].addFace(xFace);
//					}
//					else if(i == getNumberXCells())
//					{
//						xFace.setBoundaryFace(true);
//						getCellsArray()[i-1][j].addFace(xFace);
//					}
//					else
//					{
//						getCellsArray()[i][j].addFace(xFace);
//						getCellsArray()[i-1][j].addFace(xFace);
//					}
//				}
//				if(i < getNumberXCells())
//				{
//					TPFace yFace = new TPFace(getCellsX().get(i), getyStart() + j * hy, 1);
//					getFaces().add(yFace);
//					if (j == 0)
//					{
//						yFace.setBoundaryFace(true);
//						getCellsArray()[i][j].addFace(yFace);
//					} else if (j == getNumberYCells())
//					{
//						yFace.setBoundaryFace(true);
//						getCellsArray()[i][j - 1].addFace(yFace);
//					} else
//					{
//						getCellsArray()[i][j].addFace(yFace);
//						getCellsArray()[i][j - 1].addFace(yFace);
//					}
//				}
//			}
//		}
//	}
//
//	@Override
//	public void initializeSystemMatrix()
//	{
//		systemMatrix = new DoubleTensor(getShapeFunctions().size(), getShapeFunctions().size(), true);
//	}
//
//	@Override
//	public void initializeRhs()
//	{
//		rhs = new DoubleTensor(getShapeFunctions().size());
//	}
//
//	@Override
//	public DoubleTensor getRhs()
//	{
//		return rhs;
//	}
//
//	@Override
//	public DoubleTensor getSystemMatrix()
//	{
//		return systemMatrix;
//	}
//
//
//	@Override
//	public TPFESpace refine(Multimap<TPCell, TPCell> cellRefinedCellMapping,
//	                                                 Multimap<TPFace, TPFace> faceRefinedFaceMapping)
//	{
//		TPFESpace ret = new TPFESpace(getxStart(), getyStart(), getxEnd(), getyEnd(), getNumberXCells(),
//			getNumberYCells(), getPolynomialDegree());
//
//		return ret;
//	}
//
//	public int getPolynomialDegree()
//	{
//		return polynomialDegree;
//	}
//
//	@Override
//	public ArrayList<TPCell> getCells()
//	{
//		return cells;
//	}
//
//	@Override
//	public ArrayList<TPFace> getFaces()
//	{
//		return faces;
//	}
//
//	@Override
//	public ArrayList<TPShapeFunction> getShapeFunctions()
//	{
//		return shapeFunctions;
//	}
//
//	public ArrayList<Cell1D> getCellsX()
//	{
//		return cellsX;
//	}
//
//	public ArrayList<Cell1D> getCellsY()
//	{
//		return cellsY;
//	}
//
//	public TPCell[][] getCellsArray()
//	{
//		return cellsArray;
//	}
//
//	public double getxStart()
//	{
//		return xStart;
//	}
//
//	public double getyStart()
//	{
//		return yStart;
//	}
//
//	public double getxEnd()
//	{
//		return xEnd;
//	}
//
//	public double getyEnd()
//	{
//		return yEnd;
//	}
//
//	public int getNumberXCells()
//	{
//		return numberXCells;
//	}
//
//	public int getNumberYCells()
//	{
//		return numberYCells;
//	}
}
