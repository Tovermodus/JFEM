package mixed;

import basic.FEBaseTransformation;
import basic.LagrangeNodeFunctional;
import basic.ScalarFunction;
import basic.VectorFunction;
import com.google.common.collect.*;
import com.google.common.primitives.Doubles;
import com.google.common.primitives.Ints;
import linalg.*;
import tensorproduct.*;

import java.util.*;
import java.util.stream.Collectors;

public class NedelecSpace implements MixedFESpace<TPCell, TPFace, TPEdge, ContinuousTPShapeFunction,
	NedelecShapeFunction>
{
	
	List<List<Double>> coordinates1D;
	List<List<Cell1D>> cells1D;
	List<TPCell> cells;
	List<TPFace> faces;
	Set<TPEdge> edges;
	List<NedelecNodeFuncional> nodeFuncionals;
	Set<Integer> boundaryNodes;
	TreeMultimap<TPCell, MixedShapeFunction<TPCell, TPFace, TPEdge,ContinuousTPShapeFunction,
		NedelecShapeFunction>> supportOnCell;
	TreeMultimap<TPFace, MixedShapeFunction<TPCell, TPFace, TPEdge,ContinuousTPShapeFunction,
		NedelecShapeFunction>> supportOnFace;
	TreeMultimap<TPCell, NedelecNodeFuncional> functionalsOnCell;
	Map<List<Integer>, TPCell> lexicographicCellNumbers;
	Map<TPCell, FEBaseTransformation<VectorFunction, NedelecNodeFuncional,
			CoordinateVector, CoordinateMatrix, Tensor>> allTransformationMaps;
	BiMap<NedelecNodeFuncional, NedelecShapeFunction> functionalShapeFunctionMap;
	Set<MixedNedelecFunction> shapeFunctions;
	SparseMatrix systemMatrix;
	DenseVector rhs;
	final int dimension;
	volatile int cellcounter;
	volatile int facecounter;
	
	public NedelecSpace(CoordinateVector startCoordinates, CoordinateVector endCoordinates,
	                       List<Integer> cellsPerDimension, int polynomialDegree)
	{
		if (startCoordinates.getLength() != endCoordinates.getLength() || startCoordinates.getLength() != cellsPerDimension.size())
			throw new IllegalArgumentException();
		dimension = startCoordinates.getLength();
		cells1D = new ArrayList<>();
		coordinates1D = new ArrayList<>();
		supportOnCell = TreeMultimap.create();
		supportOnFace = TreeMultimap.create();
		functionalsOnCell = TreeMultimap.create();
		allTransformationMaps = new TreeMap<>();
		boundaryNodes = new TreeSet<>();
		functionalShapeFunctionMap = HashBiMap.create();
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
		assembleFaces();
		assembleEdges();
	}
	public void assembleFaces()
	{
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
	public void assembleEdges()
	{
		edges = new TreeSet<>();
		for(TPFace face:getFaces())
		{
			int[] tangentDimensions = new int[2];
			int subd = 0;
			for(int d = 0; d < 3; d ++)
				if(d != face.getFlatDimension())
					tangentDimensions[subd++] = d;
			edges.add(TPEdge.createEdgeFromFace(face, tangentDimensions[0],true));
			edges.add(TPEdge.createEdgeFromFace(face, tangentDimensions[0],false));
			edges.add(TPEdge.createEdgeFromFace(face, tangentDimensions[1],true));
			edges.add(TPEdge.createEdgeFromFace(face, tangentDimensions[1],false));
		}
	}
	
	@Override
	public void assembleFunctions(int polynomialDegree)
	{
		shapeFunctions = new TreeSet<>();
		nodeFuncionals = new ArrayList<>();
		for (TPCell c : getCells())
			for (int i = 0; i < 3 * polynomialDegree * Math.pow(polynomialDegree - 1, 2); i++)
				nodeFuncionals.add(new NedelecNodeFuncional(polynomialDegree, c, i));
		for (TPFace f : getFaces())
			for (int i = 0; i < 2 * polynomialDegree * (polynomialDegree - 1); i++)
				nodeFuncionals.add(new NedelecNodeFuncional(polynomialDegree, f, i));
		for (TPEdge e : getEdges())
			for (int i = 0; i < polynomialDegree; i++)
				nodeFuncionals.add(new NedelecNodeFuncional(polynomialDegree, e, i));
		for (NedelecNodeFuncional n : nodeFuncionals)
		{
			//System.out.println(n.getCells().size()+" "+n.getType());
			for (TPCell c : n.getCells())
				this.functionalsOnCell.put(c, n);
		}
		//for(TPCell c:getCells())
		//	System.out.println(this.functionalsOnCell.get(c).size());
		assemblePressureFunctions(polynomialDegree);
		assembleVelocityFunctions(polynomialDegree);
		int i = 0;
		for(MixedNedelecFunction shapeFunction:shapeFunctions)
			shapeFunction.setGlobalIndex(i++);
	}
	
	private Collection<TPEdge> getEdges()
	{
		return edges;
	}
	
	private void assemblePressureFunctions(int polynomialDegree)
	{
		
		for (TPCell cell : cells)
		{
			for (int i = 0; i < Math.pow(polynomialDegree , dimension); i++)
			{
				MixedNedelecFunction shapeFunction = new MixedNedelecFunction(new ContinuousTPShapeFunction(cell,
					polynomialDegree-1, i));
				shapeFunction.setGlobalIndex(shapeFunctions.size());
				shapeFunctions.add(shapeFunction);
				for(TPCell ce: shapeFunction.getCells())
					supportOnCell.put(ce, shapeFunction);
				for (TPFace face : shapeFunction.getFaces())
					supportOnFace.put(face, shapeFunction);
			}
		}
	}
	private void assembleVelocityFunctions(int polynomialDegree)
	{
		for(NedelecNodeFuncional n: nodeFuncionals)
		{
			if (n.getType() == NedelecNodeFuncional.EDGETYPE)
			{
				NedelecShapeFunction shapeFunction = new NedelecShapeFunction(n.getE(), n);
				functionalShapeFunctionMap.put(n, shapeFunction);
			}
			if (n.getType() == NedelecNodeFuncional.FACETYPE)
			{
				NedelecShapeFunction shapeFunction = new NedelecShapeFunction(n.getF(), n);
				functionalShapeFunctionMap.put(n, shapeFunction);
			}
			if (n.getType() == NedelecNodeFuncional.CELLTYPE)
			{
				NedelecShapeFunction shapeFunction = new NedelecShapeFunction(n.getC(), n);
				functionalShapeFunctionMap.put(n, shapeFunction);
			}
		}
		for (TPCell cell : cells)
		{
			FEBaseTransformation<VectorFunction, NedelecNodeFuncional,
				CoordinateVector, CoordinateMatrix, Tensor> transformationMap =
				new FEBaseTransformation<>(assembleOriginalFunctionSpaceOnCell(polynomialDegree,
					cell), new ArrayList<>(functionalsOnCell.get(cell)));
			allTransformationMaps.put(cell, transformationMap);
			for(NedelecNodeFuncional n: functionalsOnCell.get(cell))
			{
				NedelecShapeFunction shapeFunction = functionalShapeFunctionMap.get(n);
				shapeFunction.addTransformationMap(transformationMap, cell);
				MixedNedelecFunction mixedNedelecFunction = new MixedNedelecFunction(shapeFunction);
				shapeFunctions.add(mixedNedelecFunction);
				supportOnCell.put(cell, mixedNedelecFunction);
				for(TPFace f: cell.getFaces())
					supportOnFace.put(f,mixedNedelecFunction);
				
			}
			
		}
//		for(NedelecNodeFuncional n: functionalShapeFunctionMap.keySet())
//		{
//			System.out.println(n.evaluate(functionalShapeFunctionMap.get(n)));
//		}
		
	}
	private List<VectorFunction> assembleOriginalFunctionSpaceOnCell(int polynomialDegree, TPCell cell)
	{
		List<VectorFunction> functions = new ArrayList<>();
		for (int i = 0; i < Math.pow(polynomialDegree + 1, dimension-1)*(polynomialDegree) * dimension; i++)
		{
			DGNodalNedelecShapeFunction shapeFunction =new DGNodalNedelecShapeFunction(cell,
				polynomialDegree-1, i);
			functions.add(shapeFunction);
		}
		return functions;
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
	public Map<Integer, MixedShapeFunction<TPCell, TPFace, TPEdge,ContinuousTPShapeFunction,
		NedelecShapeFunction>> getShapeFunctions()
	{
		Map<Integer,
			MixedShapeFunction<TPCell,TPFace,TPEdge,ContinuousTPShapeFunction,NedelecShapeFunction>> functionNumbers = new TreeMap<>();
		for (MixedShapeFunction<TPCell,TPFace,TPEdge,ContinuousTPShapeFunction,NedelecShapeFunction> shapeFunction :
			shapeFunctions)
			functionNumbers.put(shapeFunction.getGlobalIndex(), shapeFunction);
		return functionNumbers;
	}
	
	@Override
	public Set<Integer> getFixedNodeIndices()
	{
		return boundaryNodes;
	}
	
	@Override
	public List<TPFace> getFaces()
	{
		return faces;
	}
	
	@Override
	public Collection<MixedShapeFunction<TPCell, TPFace, TPEdge,ContinuousTPShapeFunction,
		NedelecShapeFunction>> getShapeFunctionsWithSupportOnCell(TPCell cell)
	{
		return supportOnCell.get(cell);
	}
	
	@Override
	public Collection<MixedShapeFunction<TPCell, TPFace, TPEdge,ContinuousTPShapeFunction,
		NedelecShapeFunction>> getShapeFunctionsWithSupportOnFace(TPFace face)
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
	
	public void setVelocityBoundaryValues(VectorFunction boundaryValues)
	{
		MixedFunction boundaryMixed = new MixedFunction(boundaryValues);
		int progress = 0;
		for(MixedShapeFunction<TPCell,TPFace,TPEdge,ContinuousTPShapeFunction, NedelecShapeFunction> shapeFunction :
			getShapeFunctions().values())
		{
			if(!shapeFunction.isVelocity())
				continue;
			boolean boundaryShapeFunction = false;
			NedelecNodeFuncional nodeFuncional =
				shapeFunction.getVelocityShapeFunction().nodeFuncional;
			if(nodeFuncional.getType() == NedelecNodeFuncional.FACETYPE)
				if(nodeFuncional.getF().isBoundaryFace())
					boundaryShapeFunction = true;
			if(nodeFuncional.getType() == NedelecNodeFuncional.EDGETYPE)
				if(nodeFuncional.getE().isBoundaryEdge())
					boundaryShapeFunction = true;
			if(boundaryShapeFunction)
			{
				double nodeValue = shapeFunction.getNodeFunctional().evaluate(boundaryMixed);
				int shapeFunctionIndex = shapeFunction.getGlobalIndex();
				boundaryNodes.add(shapeFunctionIndex);
				systemMatrix.deleteLine(shapeFunctionIndex);
				getSystemMatrix().set(1, shapeFunctionIndex, shapeFunctionIndex);
				getRhs().set(nodeValue, shapeFunctionIndex);
			}
		}
	}
	public void setPressureBoundaryValues(ScalarFunction boundaryValues)
	{
		MixedFunction boundaryMixed = new MixedFunction(boundaryValues);
		int progress = 0;
		for(MixedShapeFunction<TPCell,TPFace,TPEdge,ContinuousTPShapeFunction, NedelecShapeFunction> shapeFunction :
			getShapeFunctions().values())
		{
			if(!shapeFunction.isPressure())
				continue;
			boolean boundaryShapeFunction = false;
			for(TPFace f: shapeFunction.getFaces())
			{
				if(f.isBoundaryFace() && f.isOnFace(((LagrangeNodeFunctional)shapeFunction.getPressureShapeFunction().getNodeFunctional()).getPoint()))
				{
					boundaryShapeFunction = true;
					break;
				}
			}
			if(boundaryShapeFunction)
			{
				System.out.println(shapeFunction.getGlobalIndex());
				double nodeValue = shapeFunction.getNodeFunctional().evaluate(boundaryMixed);
				int shapeFunctionIndex = shapeFunction.getGlobalIndex();
				boundaryNodes.add(shapeFunctionIndex);
				systemMatrix.deleteLine(shapeFunctionIndex);
				getSystemMatrix().set(1, shapeFunctionIndex, shapeFunctionIndex);
				getRhs().set(nodeValue, shapeFunctionIndex);
			}
		}
		
	}
}
