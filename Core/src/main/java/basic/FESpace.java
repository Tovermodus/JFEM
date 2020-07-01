package basic;

import com.google.common.collect.Lists;
import com.google.common.collect.Multimap;
import linalg.DoubleTensor;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ForkJoinPool;

public class FESpace
{
	protected ArrayList<Cell> cells;
	protected ArrayList<Face> faces;
	protected ArrayList<ScalarShapeFunction> shapeFunctions;
	protected DoubleTensor systemMatrix;
	protected DoubleTensor rhs;
	public FESpace()
	{
		cells = new ArrayList<>();
		faces = new ArrayList<>();
		shapeFunctions = new ArrayList<>();
	}
	
	public ArrayList<Cell> getCells()
	{
		return cells;
	}
	
	public ArrayList<ScalarShapeFunction> getShapeFunctions()
	{
		return shapeFunctions;
	}
	
	public ArrayList<Face> getFaces()
	{
		return faces;
	}
	
	public DoubleTensor getRhs()
	{
		return rhs;
	}
	
	public DoubleTensor getSystemMatrix()
	{
		return systemMatrix;
	}
	
	
	public FESpace refineAll(Multimap<Cell,Cell> cellRefinedCellMapping,
	                         Multimap<Face,Face> faceRefinedFaceMapping)
	{
		FESpace ret = new FESpace();
		ArrayList<Cell> refinedCells = ret.cells;
		ArrayList<Face> refinedFaces = ret.faces;
		ArrayList<ScalarShapeFunction> refinedShapeFunctions = ret.shapeFunctions;
		for(Cell cell:cells)
		{
			ArrayList<Cell> refCells = cell.refine(refinedFaces);
			for(Cell refCell: refCells)
			{
				refinedCells.add(refCell);
				cellRefinedCellMapping.put(cell, refCell);
			}
		}
		for(Face f:faces)
		{
			ArrayList<Face> refFaces = f.refine(cellRefinedCellMapping);
			for(Face refFace:refFaces)
			{
				refinedFaces.add(refFace);
				faceRefinedFaceMapping.put(f,refFace);
			}
		}
		for(Cell refinedCell: refinedCells)
		{
			refinedCell.distributeFunctions(refinedShapeFunctions);
		}
		return ret;

	}

	public void evaluateCellIntegrals(ArrayList<CellIntegral> cellIntegrals, ArrayList<RightHandSideIntegral> rightHandSideIntegrals)
	{
		systemMatrix = new DoubleTensor(shapeFunctions.size(),shapeFunctions.size(),true);
		rhs = new DoubleTensor(shapeFunctions.size());

		List<List<Cell>> smallerList = Lists.partition(cells,12);
		int i = 0;
		ForkJoinPool pool = new ForkJoinPool(4);
		//pool.submit(()-> dd
		smallerList.stream().parallel().forEach(smallList->
		{
			for (Cell K : smallList)
			{
				//System.out.println("evaluate cell integrals: "+(int)((1.0*i++)/(cells.size())*100)+"%");
				for (ScalarShapeFunction v : K.getShapeFunctions())
				{
					for (ScalarShapeFunction u : K.getShapeFunctions())
					{
						double integral = 0;
						for (CellIntegral cellIntegral : cellIntegrals)
						{
							integral += cellIntegral.evaluateCellIntegral(K, u, v);
						}
						if(integral != 0)
							systemMatrix.add(v.getGlobalIndex(), u.getGlobalIndex(), integral);
					}
					double integral = 0;
					for (RightHandSideIntegral rightHandSideIntegral : rightHandSideIntegrals)
					{
						integral += rightHandSideIntegral.evaluateRightHandSideIntegral(K, v);
					}
					if(integral != 0)
						rhs.add(v.getGlobalIndex(), integral);

				}
			}

		});
	}
	public void evaluateFaceIntegrals(ArrayList<FaceIntegral> faceIntegrals, ArrayList<BoundaryFaceIntegral> boundaryFaceIntegrals)
	{

		List<List<Face>> smallerList = Lists.partition(faces,12);
		int i = 0;
		ForkJoinPool pool = new ForkJoinPool(4);
		//pool.submit(()-> dd
		smallerList.stream().parallel().forEach(smallList->
		{
			for (Face F : smallList)
			{
				//System.out.println("evaluate face integrals: " + (int) ((1.0 * i++) / (faces.size()
				// ) * 100) + "%");
				for (ScalarShapeFunction u : F.getShapeFunctions())
				{
					for (ScalarShapeFunction v : F.getShapeFunctions())
					{

						double integral = 0;
						for (FaceIntegral faceIntegral : faceIntegrals)
						{
							integral += faceIntegral.evaluateFaceIntegral(F,u,v);
						}
						if(integral != 0)
							systemMatrix.add(v.getGlobalIndex(), u.getGlobalIndex(), integral);
					}
				}
				if (F.isBoundaryFace())
				{
					for (ScalarShapeFunction v : F.getShapeFunctions())
					{
						double integral = 0;
						for (BoundaryFaceIntegral boundaryFaceIntegral : boundaryFaceIntegrals)
						{
							integral +=
								boundaryFaceIntegral.evaluateBoundaryFaceIntegral(F, v);
						}
						if(integral != 0)
							rhs.add(v.getGlobalIndex(), integral);
					}
				}
			}
		});
	}

}
