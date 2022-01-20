package mixed;

import it.unimi.dsi.fastutil.ints.Int2ObjectArrayMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectLinkedOpenHashMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectMap;
import linalg.CoordinateVector;
import linalg.Vector;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.Map;

public class MixedTPFESpaceFunction<MF extends MixedShapeFunction<TPCell, TPFace, ?, ?>>
	extends MixedFESpaceFunction<MF, TPCell, TPFace>
{
	final Int2ObjectMap<HashSet<MF>> supportOnCell;
	final Int2ObjectMap<ArrayList<MF>> supportOnCellFast;
	
	public MixedTPFESpaceFunction(final MF[] functions, final double[] coefficients)
	{
		super(functions, coefficients);
		supportOnCell = new Int2ObjectArrayMap<>();
		supportOnCellFast = new Int2ObjectArrayMap<>();
		for (final MF function : functions)
		{
			for (final TPCell cell : function.getCells())
			{
				if (!supportOnCell.containsKey(cell.doneCode()))
					supportOnCell.put(cell.doneCode(), new HashSet<>());
				supportOnCell.get(cell.doneCode())
				             .add(function);
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
		supportOnCell = new Int2ObjectLinkedOpenHashMap<>();
		supportOnCellFast = new Int2ObjectLinkedOpenHashMap<>();
		for (final MF function : functions.values())
		{
			for (final TPCell cell : function.getCells())
			{
				if (!supportOnCell.containsKey(cell.doneCode()))
					supportOnCell.put(cell.doneCode(), new HashSet<>());
				supportOnCell.get(cell.doneCode())
				             .add(function);
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
		final MixedValue ret = new MixedValue(cell.getDimension());
		supportOnCellFast
			.get(cell.doneCode())
			.stream()
			.map(f -> f.valueInCell(pos, cell))
			.forEach(ret::addInPlace);
		return ret;
	}
	
	@Override
	public MixedGradient gradientInCell(final CoordinateVector pos, final TPCell cell)
	{
		return supportOnCellFast
			.get(cell.doneCode())
			.stream()
			.map(f -> f.gradientInCell(pos, cell))
			.reduce(new MixedGradient(pos.getLength()), MixedGradient::add);
	}
}
