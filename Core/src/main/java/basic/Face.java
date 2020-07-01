package basic;

import com.google.common.collect.Multimap;
import linalg.DoubleTensor;

import java.util.ArrayList;

public abstract class Face
{
	protected ArrayList<Cell> cells;
	protected ArrayList<ScalarShapeFunction> shapeFunctions;
	protected TensorFunction normal;
	private boolean isBoundaryFace;
	public Face(TensorFunction normal)
	{
		cells = new ArrayList<>();
		shapeFunctions = new ArrayList<>();
		this.normal = normal;
	}
	
	public ArrayList<Cell> getCells()
	{
		return cells;
	}
	
	public ArrayList<ScalarShapeFunction> getShapeFunctions()
	{
		return shapeFunctions;
	}
	
	public boolean isBoundaryFace()
	{
		return isBoundaryFace;
	}
	
	public TensorFunction getNormal()
	{
		return normal;
	}
	
	public abstract Cell getNormalDownstreamCell(DoubleTensor pos);
	public abstract Cell getNormalUpstreamCell(DoubleTensor pos);
	public abstract DoubleTensor center();
	public abstract boolean isOnFace(DoubleTensor pos);
	public double jumpInValue(ScalarShapeFunction func, DoubleTensor pos)
	{
		return func.valueInCell(pos,getNormalUpstreamCell(pos)) - func.valueInCell(pos,
			getNormalDownstreamCell(pos));
	}
	
	public DoubleTensor jumpInDerivative(ScalarShapeFunction func, DoubleTensor pos)
	{
		return func.derivativeInCell(pos,getNormalUpstreamCell(pos)).sub(
			func.derivativeInCell(pos, getNormalDownstreamCell(pos)));
	}
	
	public double averageInValue(ScalarShapeFunction func, DoubleTensor pos)
	{
		return 0.5*(func.valueInCell(pos,getNormalUpstreamCell(pos))+func.valueInCell(pos,
			getNormalDownstreamCell(pos)));
	}
	
	public DoubleTensor averageInDerivative(ScalarShapeFunction func, DoubleTensor pos)
	{
		return func.derivativeInCell(pos,getNormalUpstreamCell(pos)).add(
			func.derivativeInCell(pos,getNormalDownstreamCell(pos))).mul(0.5);
	}
	public DoubleTensor normalAverageInValue(ScalarShapeFunction func, DoubleTensor pos)
	{
		return normal.value(pos).mul(0.5*jumpInValue(func, pos));
	}
	public double normalAverageInDerivative(ScalarShapeFunction func, DoubleTensor pos)
	{
		return normal.value(pos).inner(jumpInDerivative(func, pos));
	}
	public abstract ArrayList<Face> refine(Multimap<Cell, Cell> cellMap);
	
	public void setBoundaryFace(boolean boundaryFace)
	{
		isBoundaryFace = boundaryFace;
	}
}
