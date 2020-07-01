package basic;

import linalg.DoubleTensor;

import java.util.ArrayList;
import java.util.Map;

public abstract class ScalarShapeFunction extends ScalarFunction
{
	protected ScalarNodeFunctional nodeFunctional;
	protected ArrayList<Cell> cells;
	private int globalIndex;
	public void setGlobalIndex(int index)
	{
		this.globalIndex = index;
	}
	
	public ArrayList<Cell> getCells()
	{
		return cells;
	}
	
	public int getGlobalIndex()
	{
		return globalIndex;
	}
	
	public ScalarNodeFunctional getNodeFunctional()
	{
		return nodeFunctional;
	}
	
	
	public abstract double valueInCell(DoubleTensor pos, Cell cell);
	public abstract DoubleTensor derivativeInCell(DoubleTensor pos, Cell cell);
	public abstract Map<Integer, Double> prolongate(ArrayList<ScalarShapeFunction> refinedFunctions);
}
