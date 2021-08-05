package tensorproduct;

import basic.Assembleable;
import basic.MatrixFESpace;
import com.google.common.collect.Iterables;
import com.google.common.collect.Lists;
import com.google.common.collect.TreeMultimap;
import com.google.common.primitives.Doubles;
import com.google.common.primitives.Ints;
import linalg.*;
import tensorproduct.geometry.Cell1D;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

import java.util.*;
import java.util.stream.Collectors;

public class ScalarRTFESpace extends CartesianGridSpace<RTComponentFunction>
{
	public ScalarRTFESpace(CoordinateVector startCoordinates, CoordinateVector endCoordinates,
	                       List<Integer> cellsPerDimension)
	{
		super(startCoordinates, endCoordinates, cellsPerDimension);
	}
	
	@Override
	public void assembleFunctions(int polynomialDegree)
	{
		shapeFunctions = new TreeSet<>();
		for(TPCell cell: getCells())
			for(int i = 0; i < Math.pow(polynomialDegree+1, getDimension()-1)*(polynomialDegree+2); i++)
			{
				RTComponentFunction shapeFunction = new RTComponentFunction(cell,
					polynomialDegree, i,0);
				shapeFunction.setGlobalIndex(shapeFunctions.size());
				shapeFunctions.add(shapeFunction);
				for(TPCell supportCell: shapeFunction.getCells())
					supportOnCell.put(supportCell, shapeFunction);
				for(TPFace supportFace: shapeFunction.getFaces())
					supportOnFace.put(supportFace, shapeFunction);
			}
		if(getDimension() == 2)
			if(shapeFunctions.size() != ((polynomialDegree+1)*grid.cellsPerDimension.get(0)+1)*(polynomialDegree+1)*grid.cellsPerDimension.get(1))
				throw new IllegalStateException("Identification did not work ");
		if(getDimension() == 3)
			if(shapeFunctions.size() != ((polynomialDegree+1)*grid.cellsPerDimension.get(0)+1)
				*(polynomialDegree+1)*grid.cellsPerDimension.get(1)
				*(polynomialDegree+1)*grid.cellsPerDimension.get(2))
				throw new IllegalStateException("Identification did not work");
	}
	
}

