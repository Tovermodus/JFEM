package tensorproduct;

import basic.VectorFESpaceFunction;
import basic.VectorFunctionOnCells;
import basic.VectorShapeFunction;
import it.unimi.dsi.fastutil.ints.Int2ObjectLinkedOpenHashMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectMap;
import linalg.CoordinateDenseMatrix;
import linalg.CoordinateMatrix;
import linalg.CoordinateVector;
import linalg.Vector;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.Map;

public class VectorTPFESpaceFunction<ST extends VectorShapeFunction<TPCell, TPFace>>
	extends VectorFESpaceFunction<ST>
	implements
	VectorFunctionOnCells<TPCell, TPFace>
{
	final Int2ObjectMap<HashSet<ST>> supportOnCell;
	final Int2ObjectMap<ArrayList<ST>> supportOnCellFast;
	
	public VectorTPFESpaceFunction(final ST[] functions, final double[] coefficients)
	{
		super(functions, coefficients);
		supportOnCell = new Int2ObjectLinkedOpenHashMap<>();
		supportOnCellFast = new Int2ObjectLinkedOpenHashMap<>();
		for (final ST function : functions)
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
	
	public VectorTPFESpaceFunction(final Map<Integer, ST> functions, final Vector coefficients)
	{
		super(functions, coefficients);
		supportOnCell = new Int2ObjectLinkedOpenHashMap<>();
		supportOnCellFast = new Int2ObjectLinkedOpenHashMap<>();
		for (final ST function : functions.values())
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
	public CoordinateVector valueInCell(final CoordinateVector pos, final TPCell cell)
	{
		return supportOnCellFast
			.get(cell.doneCode())
			.stream()
			.map(f -> f.valueInCell(pos, cell))
			.reduce(new CoordinateVector(pos.getLength()), CoordinateVector::add);
	}
	
	@Override
	public CoordinateMatrix gradientInCell(final CoordinateVector pos, final TPCell cell)
	{
		return supportOnCellFast
			.get(cell.doneCode())
			.stream()
			.map(f -> f.gradientInCell(pos, cell))
			.reduce(new CoordinateDenseMatrix(pos.getLength(), pos.getLength()), CoordinateMatrix::add);
	}
}
