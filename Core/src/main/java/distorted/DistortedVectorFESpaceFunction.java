package distorted;

import basic.VectorFESpaceFunction;
import distorted.geometry.DistortedCell;
import it.unimi.dsi.fastutil.ints.Int2ObjectLinkedOpenHashMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectMap;
import linalg.*;

import java.util.HashMap;
import java.util.Map;
import java.util.TreeSet;

public class DistortedVectorFESpaceFunction extends VectorFESpaceFunction<DistortedVectorShapeFunction> implements DistortedVectorFunction
{
	final Int2ObjectMap<TreeSet<DistortedVectorShapeFunction>> supportOnCell;
	
	public DistortedVectorFESpaceFunction(final DistortedVectorShapeFunction[] functions, final double[] coefficients)
	{
		super(functions, coefficients);
		supportOnCell = new Int2ObjectLinkedOpenHashMap<>();
		for (final DistortedVectorShapeFunction function : functions)
		{
			for (final DistortedCell cell : function.getCells())
			{
				if (!supportOnCell.containsKey(cell.doneCode())) supportOnCell.put(cell.doneCode(),
				                                                         new TreeSet<>());
				supportOnCell.get(cell.doneCode()).add(function);
			}
		}
	}
	
	public DistortedVectorFESpaceFunction(final Map<Integer, DistortedVectorShapeFunction> functions, final Vector coefficients)
	{
		super(functions, coefficients);
		supportOnCell = new Int2ObjectLinkedOpenHashMap<>();
		for (final DistortedVectorShapeFunction function : functions.values())
		{
			for (final DistortedCell cell : function.getCells())
			{
				if (!supportOnCell.containsKey(cell.doneCode())) supportOnCell.put(cell.doneCode(),
				                                                                   new TreeSet<>());
				supportOnCell.get(cell.doneCode()).add(function);
			}
		}
	}
	
	@Override
	public CoordinateVector valueOnReferenceCell(final CoordinateVector pos, final DistortedCell cell)
	{
		final CoordinateVector acc = new CoordinateVector(getDomainDimension());
		for (final DistortedVectorShapeFunction distortedVectorShapeFunction : supportOnCell.get(cell.doneCode()))
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
		for (final DistortedVectorShapeFunction distortedVectorShapeFunction : supportOnCell.get(cell.doneCode()))
		{
			final Rank1CoordinateMatrix mul = distortedVectorShapeFunction.gradientOnReferenceCell(pos,
			                                                                                       cell);
			mul.mulInPlace(coefficients.get(distortedVectorShapeFunction));
			acc.addInPlace(mul);
		}
		return acc;
	}
}
