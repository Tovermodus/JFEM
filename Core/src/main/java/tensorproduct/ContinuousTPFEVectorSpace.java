package tensorproduct;

import basic.*;
import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;
import com.google.common.collect.TreeMultimap;
import com.google.common.primitives.Doubles;
import com.google.common.primitives.Ints;
import linalg.*;
import linalg.Vector;

import java.util.*;
import java.util.stream.Collectors;

public class ContinuousTPFEVectorSpace implements MatrixFESpace<TPCell, TPFace, TPEdge,ContinuousTPVectorFunction,
	CoordinateVector,
	CoordinateMatrix,
	Tensor>,
	Assembleable
{
	List<List<Double>> coordinates1D;
	List<List<Cell1D>> cells1D;
	List<TPCell> cells;
	List<TPFace> faces;
	TreeMultimap<TPCell, ContinuousTPVectorFunction> supportOnCell;
	TreeMultimap<TPFace, ContinuousTPVectorFunction> supportOnFace;
	Map<List<Integer>, TPCell> lexicographicCellNumbers;
	TreeSet<ContinuousTPVectorFunction> shapeFunctions;
	SparseMatrix systemMatrix;
	DenseVector rhs;
	final int dimension;
	volatile int cellcounter;
	volatile int facecounter;
	
	public ContinuousTPFEVectorSpace(CoordinateVector startCoordinates, CoordinateVector endCoordinates,
	                                 List<Integer> cellsPerDimension, int polynomialDegree)
	{
		if (startCoordinates.getLength() != endCoordinates.getLength() || startCoordinates.getLength() != cellsPerDimension.size())
			throw new IllegalArgumentException();
		dimension = startCoordinates.getLength();
		cells1D = new ArrayList<>();
		coordinates1D = new ArrayList<>();
		supportOnCell = TreeMultimap.create();
		supportOnFace = TreeMultimap.create();
		QuadratureRule1D quad;
		if (polynomialDegree < 3)
			quad = QuadratureRule1D.Gauss5;
		else
			quad = QuadratureRule1D.Gauss5;
		for (int i = 0; i < startCoordinates.getLength(); i++)
		{
			List<Cell1D> cellsForDirection = new ArrayList<>();
			List<Double> coordinatesForDirection = new ArrayList<>();
			double le = (endCoordinates.at(i) - startCoordinates.at(i)) / cellsPerDimension.get(i);
			for (int j = 0; j < cellsPerDimension.get(i); j++)
			{
				coordinatesForDirection.add(startCoordinates.at(i) + le * j);
				cellsForDirection.add(new Cell1D(startCoordinates.at(i) + le * j,
					startCoordinates.at(i) + le * (j + 1), quad));
			}
			coordinatesForDirection.add(endCoordinates.at(i));
			cells1D.add(cellsForDirection);
			coordinates1D.add(coordinatesForDirection);
		}
		
	}
	
	private <T> int[] decomposeLexicographic(int index, List<List<T>> list)
	{
		int[] ret = new int[list.size()];
		for (int i = list.size() - 1; i >= 0; i--)
		{
			ret[i] = index % list.get(i).size();
			index = index / list.get(i).size();
		}
		if (index != 0)
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
	public void assembleFunctions(int polynomialDegree)
	{
		
		shapeFunctions = new TreeSet<>();
		for (TPCell cell : cells)
		{
			for (int i = 0; i < Math.pow(polynomialDegree + 1, dimension) * dimension; i++)
			{
				ContinuousTPVectorFunction shapeFunction = new ContinuousTPVectorFunction(cell, polynomialDegree, i,
					ContinuousTPShapeFunction.class);
				shapeFunction.setGlobalIndex(shapeFunctions.size());
				shapeFunctions.add(shapeFunction);
				for (TPCell supportCell : shapeFunction.getCells())
					supportOnCell.put(supportCell, shapeFunction);
				for (TPFace supportFace : shapeFunction.getFaces())
					supportOnFace.put(supportFace, shapeFunction);
			}
		}
		if (dimension == 2)
			if (shapeFunctions.size() != dimension * (polynomialDegree * cells1D.get(0).size() + 1) * (polynomialDegree * cells1D.get(1).size() + 1))
				throw new IllegalStateException("Identification did not work");
		if (dimension == 3)
			if (shapeFunctions.size() != dimension * (polynomialDegree * cells1D.get(0).size() + 1) * (polynomialDegree * cells1D.get(1).size() + 1) * (polynomialDegree * cells1D.get(2).size() + 1))
				throw new IllegalStateException("Identification did not work");
	}
	
	
	@Override
	public void initializeSystemMatrix()
	{
		systemMatrix = new SparseMatrix(shapeFunctions.size(), shapeFunctions.size());
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
	public synchronized int increaseCellCounter()
	{
		return cellcounter++;
	}
	
	@Override
	public synchronized int increaseFaceCounter()
	{
		return facecounter++;
	}
	
	@Override
	public SparseMatrix getSystemMatrix()
	{
		return systemMatrix;
	}
	
	@Override
	public int getDimension()
	{
		return dimension;
	}
	
	@Override
	public List<TPCell> getCells()
	{
		return cells;
	}
	
	@Override
	public Map<Integer, ContinuousTPVectorFunction> getShapeFunctions()
	{
		Map<Integer, ContinuousTPVectorFunction> functionNumbers = new TreeMap<>();
		for (ContinuousTPVectorFunction shapeFunction : shapeFunctions)
			functionNumbers.put(shapeFunction.getGlobalIndex(), shapeFunction);
		return functionNumbers;
	}
	
	
	@Override
	public List<TPFace> getFaces()
	{
		return faces;
	}
	
	@Override
	public Collection<ContinuousTPVectorFunction> getShapeFunctionsWithSupportOnCell(TPCell cell)
	{
		return supportOnCell.get(cell);
	}
	
	@Override
	public Collection<ContinuousTPVectorFunction> getShapeFunctionsWithSupportOnFace(TPFace face)
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
				coordinatesForDirection.add(doubles.get(0) + le * j + Math.random() * 1e-5);
			}
			coordinatesForDirection.add(Iterables.getLast(doubles));
			plotCoordinates1D.add(coordinatesForDirection);
		}
		return Lists.cartesianProduct(plotCoordinates1D).stream().map(Doubles::toArray).map(CoordinateVector::fromValues).collect(Collectors.toList());
	}
	public void setBoundaryValues(VectorFunction boundaryValues)
	{
		int progress = 0;
		for (TPFace face : getFaces())
		{
			System.out.println((int)(100*progress/getFaces().size())+"%");
			progress++;
			if (face.isBoundaryFace())
			{
				for (ContinuousTPVectorFunction shapeFunction : getShapeFunctionsWithSupportOnFace(face))
				{
					double nodeValue = shapeFunction.getNodeFunctional().evaluate(boundaryValues);
					if (nodeValue != 0 || face.isOnFace(shapeFunction.getNodeFunctionalPoint()))
					{
						int shapeFunctionIndex = shapeFunction.getGlobalIndex();
						for (TPCell cell : shapeFunction.getCells())
							for (ContinuousTPVectorFunction sameSupportFunction :
								getShapeFunctionsWithSupportOnCell(cell))
								systemMatrix.set(0, shapeFunctionIndex,
									sameSupportFunction.getGlobalIndex());
						getSystemMatrix().set(1, shapeFunctionIndex, shapeFunctionIndex);
						getRhs().set(nodeValue, shapeFunctionIndex);
					}
				}
			}
		}
	}
}

