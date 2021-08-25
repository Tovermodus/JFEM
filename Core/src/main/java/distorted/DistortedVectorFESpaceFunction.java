package distorted;

import basic.VectorFESpaceFunction;
import distorted.geometry.DistortedCell;
import linalg.CoordinateMatrix;
import linalg.CoordinateVector;
import linalg.Vector;

import java.util.HashMap;
import java.util.Map;
import java.util.TreeSet;

public class DistortedVectorFESpaceFunction extends VectorFESpaceFunction<DistortedVectorShapeFunction> implements DistortedVectorFunction
{
	final HashMap<DistortedCell, TreeSet<DistortedVectorShapeFunction>> supportOnCell;
	
	public DistortedVectorFESpaceFunction(final DistortedVectorShapeFunction[] functions, final double[] coefficients)
	{
		super(functions, coefficients);
		supportOnCell = new HashMap<>();
		for (final DistortedVectorShapeFunction function : functions)
		{
			for (final DistortedCell cell : function.getCells())
			{
				if (!supportOnCell.containsKey(cell)) supportOnCell.put(cell, new TreeSet<>());
				supportOnCell.get(cell).add(function);
			}
		}
	}
	
	public DistortedVectorFESpaceFunction(final Map<Integer, DistortedVectorShapeFunction> functions, final Vector coefficients)
	{
		super(functions, coefficients);
		supportOnCell = new HashMap<>();
		for (final DistortedVectorShapeFunction function : functions.values())
		{
			for (final DistortedCell cell : function.getCells())
			{
				if (!supportOnCell.containsKey(cell)) supportOnCell.put(cell, new TreeSet<>());
				supportOnCell.get(cell).add(function);
			}
		}
	}
	
	public CoordinateVector valueOnReferenceCell(final CoordinateVector pos, final DistortedCell cell)
	{
		return supportOnCell
			.get(cell)
			.stream()
			.map(distortedVectorShapeFunction -> distortedVectorShapeFunction
				.valueOnReferenceCell(pos, cell)
				.mul(coefficients.get(distortedVectorShapeFunction)))
			.reduce(new CoordinateVector(getDomainDimension()), CoordinateVector::add);
	}
	
	public CoordinateMatrix gradientOnReferenceCell(final CoordinateVector pos, final DistortedCell cell)
	{
		return supportOnCell
			.get(cell)
			.stream()
			.map(distortedVectorShapeFunction -> distortedVectorShapeFunction
				.gradientOnReferenceCell(pos, cell)
				.mul(coefficients.get(distortedVectorShapeFunction)))
			.reduce(new CoordinateMatrix(getRangeDimension(), getDomainDimension()), CoordinateMatrix::add);
	}
}
