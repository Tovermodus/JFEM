package mixed;

import basic.FEBaseTransformation;
import basic.VectorFunction;
import basic.VectorShapeFunction;
import linalg.CoordinateMatrix;
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
public class NedelecShapeFunction extends VectorShapeFunction<TPCell, TPFace, TPEdge, NedelecShapeFunction>
{
	Set<TPCell> cells;
	Set<TPFace> faces;
	NedelecNodeFuncional nodeFuncional;
	Map<TPCell, FEBaseTransformation<VectorFunction, NedelecNodeFuncional, CoordinateVector, CoordinateMatrix, Tensor>> transformationMap;
	
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
		CoordinateVector, CoordinateMatrix, Tensor> transformation, TPCell cell)
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
	
	@Override
	public void addFace(TPFace face)
	{
		faces.add(face);
	}
	
	@Override
	public void addCell(TPCell cell)
	{
		cells.add(cell);
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
	
	@Override
	public int compareTo(@NotNull NedelecShapeFunction o)
	{
		return nodeFuncional.compareTo(o.nodeFuncional);
	}
}
