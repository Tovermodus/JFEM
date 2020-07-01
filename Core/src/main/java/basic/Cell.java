package basic;

import linalg.DoubleTensor;

import java.util.ArrayList;

public abstract class Cell
{
	protected ArrayList<ScalarShapeFunction> shapeFunctions;
	protected ArrayList<Face> faces;
	private boolean refined;
	public Cell()
	{
		faces = new ArrayList<>();
		shapeFunctions = new ArrayList<>();
	}
	public ArrayList<ScalarShapeFunction> getShapeFunctions()
	{
		return shapeFunctions;
	}
	
	public ArrayList<Face> getFaces()
	{
		return faces;
	}
	
	public boolean isRefined()
	{
		return refined;
	}
	
	public void setRefined(boolean refined)
	{
		this.refined = refined;
	}
	
	public abstract void addFace(Face face);
	public abstract void distributeFunctions(ArrayList<ScalarShapeFunction> shapeFunctions);
	public abstract boolean isInCell(DoubleTensor pos);
	public abstract DoubleTensor center();
	public abstract ArrayList<Cell> refine(ArrayList<Face> refinedFaces);

}
