package distorted;

import basic.VectorFESpaceFunction;
import distorted.geometry.DistortedCell;
import linalg.*;

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
	
	@Override
	public CoordinateVector valueOnReferenceCell(final CoordinateVector pos, final DistortedCell cell)
	{
		final CoordinateVector acc = new CoordinateVector(getDomainDimension());
		for (final DistortedVectorShapeFunction distortedVectorShapeFunction : supportOnCell.get(cell))
		{
			final CoordinateVector mul = distortedVectorShapeFunction.valueOnReferenceCell(pos, cell);
			mul.mulInPlace(coefficients.get(distortedVectorShapeFunction));
			acc.addInPlace(mul);
		}
		return acc;
	}
	
	@Override
	public CoordinateMatrix gradientOnReferenceCell(final CoordinateVector pos, final DistortedCell cell)
	{
		final CoordinateDenseMatrix acc = new CoordinateDenseMatrix(getRangeDimension(), getDomainDimension());
		for (final DistortedVectorShapeFunction distortedVectorShapeFunction : supportOnCell.get(cell))
		{
			final Rank1CoordinateMatrix mul = distortedVectorShapeFunction.gradientOnReferenceCell(pos,
			                                                                                       cell);
			mul.mulInPlace(coefficients.get(distortedVectorShapeFunction));
			acc.addInPlace(mul);
		}
		return acc;
	}
}
