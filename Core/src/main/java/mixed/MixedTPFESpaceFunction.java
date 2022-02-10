package mixed;

import io.vavr.Tuple2;
import it.unimi.dsi.fastutil.ints.Int2ObjectArrayMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectMap;
import linalg.CoordinateVector;
import linalg.Vector;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

import java.util.ArrayList;
import java.util.Map;
import java.util.TreeSet;

public class MixedTPFESpaceFunction<MF extends MixedShapeFunction<TPCell, TPFace, ?, ?>>
	extends MixedFESpaceFunction<MF, TPCell, TPFace>
{
	Int2ObjectMap<TreeSet<Tuple2<MF, Double>>> supportOnCell;
	Int2ObjectMap<ArrayList<Tuple2<MF, Double>>> supportOnCellFast;
	
	public MixedTPFESpaceFunction(final MF[] functions, final double[] coefficients)
	{
		super(functions, coefficients);
		supportOnCell = new Int2ObjectArrayMap<>();
		supportOnCellFast = new Int2ObjectArrayMap<>();
		for (final var function : functions)
		{
			for (final TPCell cell : function.getCells())
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
	
	public MixedTPFESpaceFunction(final Map<Integer, MF> functions, final Vector coefficients)
	{
		super(functions, coefficients);
		supportOnCell = new Int2ObjectArrayMap<>();
		supportOnCellFast = new Int2ObjectArrayMap<>();
		for (final var function : functions.values())
		{
			for (final TPCell cell : function.getCells())
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
	public MixedValue valueInCell(final CoordinateVector pos, final TPCell cell)
	{
		final MixedValue acc = new MixedValue(cell.getDimension());
		for (final var shapeFunction :
			supportOnCellFast.get(
				cell.doneCode()))
		{
			final MixedValue mul = shapeFunction._1.valueInCell(pos, cell);
			mul.mulInPlace(shapeFunction._2);
			acc.addInPlace(mul);
		}
		return acc;
	}
	
	@Override
	public MixedGradient gradientInCell(final CoordinateVector pos, final TPCell cell)
	{
		final MixedGradient acc = new MixedGradient(cell.getDimension());
		for (final var shapeFunction :
			supportOnCellFast.get(
				cell.doneCode()))
		{
			final MixedGradient mul = shapeFunction._1.gradientInCell(pos, cell);
			mul.mulInPlace(shapeFunction._2);
			acc.addInPlace(mul);
		}
		return acc;
	}
}
