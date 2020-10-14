package tensorproduct;

import basic.FEBaseTransformation;
import basic.NodeFunctional;
import basic.VectorFunction;
import basic.VectorShapeFunction;
import linalg.CoordinateMatrix;
import linalg.CoordinateVector;
import linalg.Tensor;
import mixed.NedelecNodeFuncional;
import org.jetbrains.annotations.NotNull;

import java.util.Map;
import java.util.Set;
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
		nodeFuncional = funcional;
	}
	public NedelecShapeFunction(TPFace face, NedelecNodeFuncional funcional)
	{
		cells = new TreeSet<>();
		faces = new TreeSet<>();
		faces.add(face);
		cells.addAll(face.getCells());
		nodeFuncional = funcional;
		
	}
	public NedelecShapeFunction(TPEdge edge, NedelecNodeFuncional funcional)
	{
		cells = new TreeSet<>();
		faces = new TreeSet<>();
		cells.addAll(edge.getCells());
		faces.addAll(edge.getFaces());
		nodeFuncional = funcional;
		
		
	}
	public void addTransformationMaps(Map<TPCell, FEBaseTransformation<VectorFunction, NedelecNodeFuncional,
		CoordinateVector, CoordinateMatrix, Tensor>> allTransformationMaps)
	{
		for(TPCell cell:getCells())
			transformationMap.put(cell,allTransformationMaps.get(cell));
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
		if(nodeFuncional.getType()>o.nodeFuncional.getType())
			return 1;
		if(nodeFuncional.getType()<o.nodeFuncional.getType())
			return -1;
		if(nodeFuncional.getType() == NedelecNodeFuncional.EDGETYPE)
			return nodeFuncional.getE().compareTo(o.nodeFuncional.getE());
		if(nodeFuncional.getType() == NedelecNodeFuncional.FACETYPE)
			return nodeFuncional.getF().compareTo(o.nodeFuncional.getF());
		if(nodeFuncional.getType() == NedelecNodeFuncional.CELLTYPE)
			return nodeFuncional.getC().compareTo(o.nodeFuncional.getC());
		return  1;
	}
}
