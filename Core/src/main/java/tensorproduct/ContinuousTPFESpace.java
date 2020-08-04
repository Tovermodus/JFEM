package tensorproduct;

import basic.Assembleable;
import basic.MatrixFESpace;
import basic.VectorFunction;
import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;
import com.google.common.primitives.Doubles;
import com.google.common.primitives.Ints;
import linalg.*;
import linalg.Vector;

import java.util.*;
import java.util.stream.Collectors;

public class ContinuousTPFESpace implements MatrixFESpace<TPCell<ContinuousTPShapeFunction>,
	TPFace<ContinuousTPShapeFunction>,ContinuousTPShapeFunction,Double, Vector, Matrix,
	ContinuousTPFESpace>, Assembleable
{
	List<List<Double>> coordinates1D;
	List<List<Cell1D>> cells1D;
	List<TPCell<ContinuousTPShapeFunction>> cells;
	List<TPFace<ContinuousTPShapeFunction>> faces;
	Map<List<Integer>, TPCell<ContinuousTPShapeFunction>> lexicographicCellNumbers;
	Set<ContinuousTPShapeFunction> shapeFunctions;
	SparseMatrix systemMatrix;
	DenseVector rhs;
	int dimension;
	
	public ContinuousTPFESpace(CoordinateVector startCoordinates, CoordinateVector endCoordinates,
	                 List<Integer> cellsPerDimension, int polynomialDegree)
	{
		if (startCoordinates.getLength() != endCoordinates.getLength() || startCoordinates.getLength() != cellsPerDimension.size())
			throw new IllegalArgumentException();
		cells1D = new ArrayList<>();
		coordinates1D = new ArrayList<>();
		dimension = startCoordinates.getLength();
		QuadratureRule1D quad;
		if (polynomialDegree < 3)
			quad = QuadratureRule1D.Gauss5;
		else
			quad = QuadratureRule1D.Gauss5;
		for (int i = 0; i < dimension; i++)
		{
			List<Cell1D> cellsForDirection = new ArrayList<>();
			List<Double> coordinatesForDirection = new ArrayList<>();
			double le = (endCoordinates.at(i) - startCoordinates.at(i)) / cellsPerDimension.get(i);
			for (int j = 0; j < cellsPerDimension.get(i); j++)
			{
				coordinatesForDirection.add(startCoordinates.at(i) + le * j);
				cellsForDirection.add(new Cell1D(startCoordinates.at(i) + le * j,
					startCoordinates.at(i) + le * (j + 1), quad, cellsForDirection.size()));
			}
			coordinatesForDirection.add(endCoordinates.at(i));
			cells1D.add(cellsForDirection);
			coordinates1D.add(coordinatesForDirection);
		}
		
	}
	
	@Override
	public void assembleFunctions(int polynomialDegree)
	{
		shapeFunctions = new TreeSet<>();
		for(TPCell<ContinuousTPShapeFunction> cell: cells)
			for(int i = 0; i < Math.pow(polynomialDegree+1, dimension); i++)
			{
				ContinuousTPShapeFunction shapeFunction = new ContinuousTPShapeFunction(cell, i,
					polynomialDegree);
				shapeFunction.setGlobalIndex(shapeFunctions.size());
				shapeFunctions.add(shapeFunction);
				
			}
		if(dimension == 2)
			if(shapeFunctions.size() != (polynomialDegree*cells1D.get(0).size()+1)*(polynomialDegree*cells1D.get(1).size()+1))
				throw new IllegalStateException("Identification did not work");
		if(dimension == 3)
			if(shapeFunctions.size() != (polynomialDegree*cells1D.get(0).size()+1)*(polynomialDegree*cells1D.get(1).size()+1)*(polynomialDegree*cells1D.get(2).size()+1))
				throw new IllegalStateException("Identification did not work");
	}
	
	@Override
	public void assembleCells()
	{
		cells = new ArrayList<>();
		faces = new ArrayList<>();
		lexicographicCellNumbers = new HashMap<>();
		for (List<Cell1D> cellCombination : Lists.cartesianProduct(cells1D))
		{
			cells.add(new TPCell<>(cellCombination));
			lexicographicCellNumbers.put(cellCombination.stream().map(Cell1D::getIndexInDimension).collect(Collectors.toList()),
				Iterables.getLast(cells));
		}
		for (int i = 0; i < dimension; i++)
		{
			List<List<Cell1D>> cells1Dcopy = new ArrayList<>(cells1D);
			cells1Dcopy.remove(i);
			for (int j = 0; j < coordinates1D.get(i).size(); j++)
			{
				for (List<Cell1D> cellCombination : Lists.cartesianProduct(cells1Dcopy))
				{
					faces.add(new TPFace<>(cellCombination, i, coordinates1D.get(i).get(j),
						(j == 0 || j == coordinates1D.get(i).size() - 1)));
					int[] lexicographicIndex = new int[coordinates1D.size()];
					int[] subLexIndex =
						cellCombination.stream().mapToInt(Cell1D::getIndexInDimension).toArray();
					int subd = 0;
					for (int k = 0; k < coordinates1D.size(); k++)
					{
						if (k == i)
							lexicographicIndex[k] = j - 1;
						else
							lexicographicIndex[k] = subLexIndex[subd++];
					}
					if (j != 0)
					{
						lexicographicCellNumbers.get(Ints.asList(lexicographicIndex)).addFace(Iterables.getLast(faces));
					}
					if (j != cells1D.get(i).size())
					{
						lexicographicIndex[i] = j;
						lexicographicCellNumbers.get(Ints.asList(lexicographicIndex)).addFace(Iterables.getLast(faces));
					}
					
				}
			}
		}
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
	public SparseMatrix getSystemMatrix()
	{
		return systemMatrix;
	}
	
	@Override
	public List<TPCell<ContinuousTPShapeFunction>> getCells()
	{
		return cells;
	}
	
	@Override
	public Map<Integer, ContinuousTPShapeFunction> getShapeFunctions()
	{
		Map<Integer, ContinuousTPShapeFunction> functionNumbers = new TreeMap<>();
		for(ContinuousTPShapeFunction shapeFunction:shapeFunctions)
			functionNumbers.put(shapeFunction.getGlobalIndex(), shapeFunction);
		return functionNumbers;
	}
	
	@Override
	public List<TPFace<ContinuousTPShapeFunction>> getFaces()
	{
		return faces;
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
	
}

