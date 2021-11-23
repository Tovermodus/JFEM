package distorted;

import basic.VectorFESpaceFunction;
import distorted.geometry.DistortedCell;
import io.vavr.Tuple2;
import it.unimi.dsi.fastutil.ints.Int2ObjectLinkedOpenHashMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectMap;
import linalg.*;

import java.util.ArrayList;
import java.util.Map;
import java.util.TreeSet;

public class DistortedVectorFESpaceFunction
	extends VectorFESpaceFunction<DistortedVectorShapeFunction>
	implements DistortedVectorFunctionOnCells
{
	Int2ObjectMap<TreeSet<Tuple2<DistortedVectorShapeFunction, Double>>> supportOnCell;
	Int2ObjectMap<ArrayList<Tuple2<DistortedVectorShapeFunction, Double>>> supportOnCellFast;
	
	public DistortedVectorFESpaceFunction(final DistortedVectorShapeFunction[] functions,
	                                      final double[] coefficients)
	{
		super(functions, coefficients);
		supportOnCell = new Int2ObjectLinkedOpenHashMap<>();
		supportOnCellFast = new Int2ObjectLinkedOpenHashMap<>();
		for (final DistortedVectorShapeFunction function : functions)
		{
			for (final DistortedCell cell : function.getCells())
			{
				if (!supportOnCell.containsKey(cell.doneCode()))
					supportOnCell.put(cell.doneCode(), new TreeSet<>());
				supportOnCell.get(cell.doneCode())
				             .add(new Tuple2<>(function, coefficients[function.getGlobalIndex()]));
			}
		}
		for (final int cellCode : supportOnCell.keySet())
		{
			supportOnCellFast.put(cellCode, new ArrayList<>(supportOnCell.get(cellCode)));
		}
	}
	
	public DistortedVectorFESpaceFunction(final Map<Integer, DistortedVectorShapeFunction> functions,
	                                      final Vector coefficients)
	{
		super(functions, coefficients);
		resetCoefficients(functions, coefficients);
	}
	
	@Override
	public void resetCoefficients(final Map<Integer, DistortedVectorShapeFunction> functions,
	                              final Vector coefficients)
	{
		supportOnCell = new Int2ObjectLinkedOpenHashMap<>();
		supportOnCellFast = new Int2ObjectLinkedOpenHashMap<>();
		for (final DistortedVectorShapeFunction function : functions.values())
		{
			for (final DistortedCell cell : function.getCells())
			{
				if (!supportOnCell.containsKey(cell.doneCode()))
					supportOnCell.put(cell.doneCode(), new TreeSet<>());
				supportOnCell.get(cell.doneCode())
				             .add(new Tuple2<>(function, coefficients.at(function.getGlobalIndex())));
			}
		}
		for (final int cellCode : supportOnCell.keySet())
		{
			supportOnCellFast.put(cellCode, new ArrayList<>(supportOnCell.get(cellCode)));
		}
	}
	
	@Override
	public CoordinateVector valueOnReferenceCell(final CoordinateVector pos, final DistortedCell cell)
	{
		final CoordinateVector acc = new CoordinateVector(getDomainDimension());
		for (final Tuple2<DistortedVectorShapeFunction, Double> distortedVectorShapeFunction :
			supportOnCellFast.get(
				cell.doneCode()))
		{
			final CoordinateVector mul = distortedVectorShapeFunction._1.valueOnReferenceCell(pos, cell);
			mul.mulInPlace(distortedVectorShapeFunction._2);
			acc.addInPlace(mul);
		}
		return acc;
	}
	
	@Override
	public CoordinateMatrix gradientOnReferenceCell(final CoordinateVector pos, final DistortedCell cell)
	{
		final CoordinateDenseMatrix acc = new CoordinateDenseMatrix(getRangeDimension(), getDomainDimension());
		for (final Tuple2<DistortedVectorShapeFunction, Double> distortedVectorShapeFunction : supportOnCellFast.get(
			cell.doneCode()))
		{
			final Rank1CoordinateMatrix mul = distortedVectorShapeFunction._1.gradientOnReferenceCell(pos,
			                                                                                          cell);
			mul.mulInPlace(distortedVectorShapeFunction._2);
			acc.addInPlace(mul);
		}
		return acc;
	}
}
