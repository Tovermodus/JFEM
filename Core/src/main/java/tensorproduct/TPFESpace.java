package tensorproduct;

import basic.*;
import com.google.common.collect.*;
import com.google.common.primitives.Doubles;
import com.google.common.primitives.Ints;
import linalg.*;
import linalg.Vector;
import org.checkerframework.checker.nullness.qual.Nullable;

import java.awt.image.AreaAveragingScaleFilter;
import java.util.*;
import java.util.stream.Collectors;

/*
Workflow of this finite element is: first generate 1D cells, then
 */
public class TPFESpace implements MatrixFESpace<TPCell,
	TPFace,
	TPShapeFunction,Double, CoordinateVector,Matrix,TPFESpace>,
	Assembleable
{
	List<List<Double>> coordinates1D;
	List<List<LagrangeBasisFunction1D>> basisFunctions;
	List<List<Cell1D>> cells1D;
	List<TPCell> cells;
	List<TPFace> faces;
	Map<List<Integer>, TPCell> lexicographicCellNumbers;
	TreeMultimap<TPCell,TPShapeFunction> supportOnCell;
	TreeMultimap<TPFace,TPShapeFunction> supportOnFace;
	List<TPShapeFunction> shapeFunctions;
	SparseMatrix systemMatrix;
	DenseVector rhs;
	public TPFESpace(CoordinateVector startCoordinates, CoordinateVector endCoordinates,
	                 List<Integer> cellsPerDimension, int polynomialDegree)
	{
		if(startCoordinates.getLength() != endCoordinates.getLength()|| startCoordinates.getLength() != cellsPerDimension.size())
			throw new IllegalArgumentException();
		cells1D = new ArrayList<>();
		coordinates1D = new ArrayList<>();
		supportOnCell = TreeMultimap.create();
		supportOnFace = TreeMultimap.create();
		QuadratureRule1D quad;
		if(polynomialDegree < 3)
			quad = QuadratureRule1D.Gauss5;
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
		for (int d = 0; d < getDimension(); d++)
		{
			List<LagrangeBasisFunction1D> functionsForDirection = new ArrayList<>();
			for (Cell1D oneDimensionalCell : cells1D.get(d))
			{
				for (int deg = 0; deg < polynomialDegree + 1; deg++)
				{
					functionsForDirection.add(new LagrangeBasisFunction1D(polynomialDegree, deg,
						oneDimensionalCell));
				}
			}
			basisFunctions.add(functionsForDirection);
		}
		List<List<LagrangeBasisFunction1D>> oneDimensionalFunctionCombinations =
			Lists.cartesianProduct(basisFunctions);
		for (int i = 0; i < oneDimensionalFunctionCombinations.size(); i++)
		{
			int[] lexicographicFunctionIndex = decomposeLexicographic(i, basisFunctions);
			int[] lexicographicCellIndex = new int[getDimension()];
			for(int j = 0; j < getDimension(); j++)
				lexicographicCellIndex[j] = lexicographicFunctionIndex[j]/(polynomialDegree+1);
			TPCell supportCell = lexicographicCellNumbers.get(Ints.asList(lexicographicCellIndex));
			TPShapeFunction function = new TPShapeFunction(supportCell,
				oneDimensionalFunctionCombinations.get(i));
			shapeFunctions.add(function);
			supportOnCell.put(supportCell, function);
			for (TPFace supportFace : function.getFaces())
				supportOnFace.put(supportFace, function);
			Iterables.getLast(shapeFunctions).setGlobalIndex(shapeFunctions.size() - 1);
		}
	}
	private<T> int[] decomposeLexicographic(int index, List<List<T>> list)
	{
		int[] ret = new int[list.size()];
		for (int i = list.size()-1; i >= 0; i--)
		{
			ret[i] = index % list.get(i).size();
			index = index/list.get(i).size();
		}
		if(index != 0)
			throw new IllegalStateException("index too high");
		return ret;
	}
	@Override
	public void assembleCells()
	{
		cells = new ArrayList<>();
		faces = new ArrayList<>();
		lexicographicCellNumbers = new HashMap<>();
		List<List<Cell1D>> cellCombinations = Lists.cartesianProduct(cells1D);
		for (int combinationIndex = 0; combinationIndex < cellCombinations.size(); combinationIndex++)
		{
			cells.add(new TPCell(cellCombinations.get(combinationIndex)));
			lexicographicCellNumbers.put(Ints.asList(decomposeLexicographic(combinationIndex, cells1D)),
				Iterables.getLast(cells));
		}
		for (int d = 0; d < getDimension(); d++)
		{
			List<List<Cell1D>> cells1Dcopy = new ArrayList<>(cells1D);
			//noinspection SuspiciousListRemoveInLoop
			cells1Dcopy.remove(d);
			for (int coordinateIndex = 0; coordinateIndex < coordinates1D.get(d).size(); coordinateIndex++)
			{
				List<List<Cell1D>> cellSubCombinations = Lists.cartesianProduct(cells1Dcopy);
				for (int combinationIndex = 0; combinationIndex < cellSubCombinations.size(); combinationIndex++)
				{
					faces.add(new TPFace(cellSubCombinations.get(combinationIndex), d,
						coordinates1D.get(d).get(coordinateIndex),
						(coordinateIndex == 0 || coordinateIndex == coordinates1D.get(d).size() - 1)));
					int[] lexicographicIndex = new int[coordinates1D.size()];
					int[] subLexIndex = decomposeLexicographic(combinationIndex, cells1Dcopy);
					int subd = 0;
					for (int k = 0; k < getDimension(); k++)
					{
						if (k == d)
							lexicographicIndex[k] = coordinateIndex - 1;
						else
							lexicographicIndex[k] = subLexIndex[subd++];
					}
					if (coordinateIndex != 0)
					{
						lexicographicCellNumbers.get(Ints.asList(lexicographicIndex)).addFace(Iterables.getLast(faces));
					}
					if (coordinateIndex != cells1D.get(d).size())
					{
						lexicographicIndex[d] = coordinateIndex;
						lexicographicCellNumbers.get(Ints.asList(lexicographicIndex)).addFace(Iterables.getLast(faces));
					}
					
				}
			}
		}
	}
	
	@Override
	public void initializeSystemMatrix()
	{
		systemMatrix = new SparseMatrix(shapeFunctions.size(),shapeFunctions.size());
	}
	
	@Override
	public void initializeRhs()
	{
		rhs = new DenseVector(shapeFunctions.size());
	}
	
	@Override
	public DenseVector getRhs()
	{
		return rhs;
	}
	
	@Override
	public SparseMatrix getSystemMatrix()
	{
		return systemMatrix;
	}
	
	@Override
	public int getDimension()
	{
		return cells1D.size();
	}
	
	@Override
	public List<TPCell> getCells()
	{
		return cells;
	}
	
	@Override
	public Map<Integer, TPShapeFunction> getShapeFunctions()
	{
		
		Map<Integer, TPShapeFunction> functionNumbers = new TreeMap<>();
		for(TPShapeFunction shapeFunction:shapeFunctions)
			functionNumbers.put(shapeFunction.getGlobalIndex(), shapeFunction);
		return functionNumbers;
	}
	
	@Override
	public List<TPFace> getFaces()
	{
		return faces;
	}
	
	@Override
	public Set<TPShapeFunction> getShapeFunctionsWithSupportOnCell(TPCell cell)
	{
		return supportOnCell.get(cell);
	}
	
	@Override
	public Set<TPShapeFunction> getShapeFunctionsWithSupportOnFace(TPFace face)
	{
		return supportOnFace.get(face);
	}
	
	@Override
	public List<CoordinateVector> generatePlotPoints(int resolution)
	{
		int n = resolution;
		List<List<Double>> plotCoordinates1D = new ArrayList<>();
		for (List<Double> doubles : coordinates1D)
		{
			List<Double> coordinatesForDirection = new ArrayList<>();
			double le =
				(Iterables.getLast(doubles) - doubles.get(0)) / n;
			for (int j = 0; j < n; j++)
			{
				coordinatesForDirection.add(doubles.get(0) + le * j+Math.random()*1e-5);
			}
			coordinatesForDirection.add(Iterables.getLast(doubles));
			plotCoordinates1D.add(coordinatesForDirection);
		}
		return Lists.cartesianProduct(plotCoordinates1D).stream().map(Doubles::toArray).map(CoordinateVector::fromValues).collect(Collectors.toList());
	}
	
	@Override
	public TPFESpace refine(Multimap<TPCell,TPCell> cellRefinedCellMapping,
	                        Multimap<TPFace,TPFace> faceRefinedFaceMapping)
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
