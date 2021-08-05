package tensorproduct;

import linalg.CoordinateVector;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

import java.util.HashSet;
import java.util.List;

public class TPVectorFESpace extends CartesianGridSpace<TPVectorFunction>
{
	public TPVectorFESpace(CoordinateVector startCoordinates, CoordinateVector endCoordinates,
	                       List<Integer> cellsPerDimension)
	{
		super(startCoordinates, endCoordinates, cellsPerDimension);
	}
	
	@Override
	public void assembleFunctions(int polynomialDegree)
	{
		shapeFunctions = new HashSet<>();
		for(TPCell cell:getCells())
		{
			for(int i = 0; i < Math.pow(polynomialDegree+1,getDimension())*getDimension(); i++)
			{
				TPVectorFunction shapeFunction = new TPVectorFunction(cell, polynomialDegree, i,
					TPShapeFunction.class);
				shapeFunction.setGlobalIndex(shapeFunctions.size());
				shapeFunctions.add(shapeFunction);
				supportOnCell.put(cell, shapeFunction);
				for(TPFace face:cell.getFaces())
					supportOnFace.put(face, shapeFunction);
			}
		}
	}
	
	
}
