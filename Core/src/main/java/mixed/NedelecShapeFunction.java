package mixed;

import basic.FEBaseTransformation;
import basic.ShapeFunction;
import basic.VectorFunction;
import basic.VectorShapeFunction;
import linalg.CoordinateMatrix;
import linalg.CoordinateTensor;
import linalg.CoordinateVector;
import linalg.Tensor;
import org.jetbrains.annotations.NotNull;
import tensorproduct.TPCell;
import tensorproduct.TPEdge;
import tensorproduct.TPFace;

import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;
public class NedelecShapeFunction implements VectorShapeFunction<TPCell, TPFace, TPEdge>, Comparable<NedelecShapeFunction>
{
	Set<TPCell> cells;
	Set<TPFace> faces;
	NedelecNodeFuncional nodeFuncional;
	Map<TPCell,
		FEBaseTransformation<VectorFunction, NedelecNodeFuncional, CoordinateVector, CoordinateMatrix, CoordinateTensor>> transformationMap;
	private int globalIndex;
	
	public NedelecShapeFunction(TPCell cell, NedelecNodeFuncional funcional)
	{
		cells = new TreeSet<>();
		faces = new TreeSet<>();
		cells.add(cell);
		faces.addAll(cell.getFaces());
		transformationMap = new TreeMap<>();
		nodeFuncional = funcional;
	}
	public NedelecShapeFunction(TPFace face, NedelecNodeFuncional funcional)
	{
		cells = new TreeSet<>();
		faces = new TreeSet<>();
		faces.add(face);
		cells.addAll(face.getCells());
		transformationMap = new TreeMap<>();
		nodeFuncional = funcional;
	}
	public NedelecShapeFunction(TPEdge edge, NedelecNodeFuncional funcional)
	{
		cells = new TreeSet<>();
		faces = new TreeSet<>();
		cells.addAll(edge.getCells());
		faces.addAll(edge.getFaces());
		transformationMap = new TreeMap<>();
		nodeFuncional = funcional;
	}
	public void addTransformationMap(FEBaseTransformation<VectorFunction, NedelecNodeFuncional,
		CoordinateVector, CoordinateMatrix, CoordinateTensor> transformation, TPCell cell)
	{
		transformationMap.put(cell,transformation);
	}
	@Override
	public Set<TPCell> getCells()
	{
		return cells;
	}
	
	@Override
	public Set<TPFace> getFaces()
	{
		return faces;
	}
	
	@Override
	public NedelecNodeFuncional getNodeFunctional()
	{
		return nodeFuncional;
	}
	
	
	public void setGlobalIndex(int index)
	{
		globalIndex = index;
	}
	@Override
	public int getGlobalIndex()
	{
		return globalIndex;
	}
	@Override
	public double divergenceInCell(CoordinateVector pos,TPCell cell)
	{
		return transformationMap.get(cell).vectorBasisFunctionDivergence(nodeFuncional,pos);
	}
	
	@Override
	public int getRangeDimension()
	{
		return getDomainDimension();
	}
	
	@Override
	public double divergence(CoordinateVector pos)
	{
		for(TPCell cell:getCells())
			if(cell.isInCell(pos))
				return divergenceInCell(pos,cell);
		return 0;
	}
	@Override
	public CoordinateVector valueInCell(CoordinateVector pos, TPCell cell)
	{
		if(!cell.isInCell(pos))
			return new CoordinateVector(pos.getLength());
		return transformationMap.get(cell).vectorBasisFunctionValue(nodeFuncional, pos);
	}
	
	@Override
	public CoordinateMatrix gradientInCell(CoordinateVector pos, TPCell cell)
	{
		if(!cell.isInCell(pos))
			return new CoordinateMatrix(pos.getLength(), pos.getLength());
		return transformationMap.get(cell).vectorBasisFunctionGradient(nodeFuncional, pos);
	}
	
	
	public int compareTo(@NotNull NedelecShapeFunction o)
	{
		return nodeFuncional.compareTo(o.nodeFuncional);
	}
	
}
