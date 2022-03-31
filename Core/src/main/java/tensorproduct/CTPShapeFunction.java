package tensorproduct;

import basic.LagrangeNodeFunctional;
import basic.PerformanceArguments;
import basic.ScalarShapeFunction;
import distorted.ReferenceCellFunction;
import linalg.CoordinateMatrix;
import linalg.CoordinateVector;
import org.jetbrains.annotations.NotNull;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

import java.util.*;
import java.util.stream.Collectors;

public class CTPShapeFunction
	implements ReferenceCellFunction<TPCell, TPFace, Double, CoordinateVector, CoordinateMatrix>,
	ScalarShapeFunction<TPCell,
		TPFace>,
	Comparable<CTPShapeFunction>

{
	private final Map<TPCell, TPShapeFunction> shapeFunctionsOnCell;
	private final int[] cellCodes;
	private final TPShapeFunction[] shapeFunctionsOnCellFast;
	private final Set<TPFace> faces;
	private final TPShapeFunction firstDefinedFunction;
	private final TPCell firstDefinedCell;
	LagrangeNodeFunctional functional;
	int polynomialDegree;
	private int globalIndex;
	
	public CTPShapeFunction(final TPCell cell, final int polynomialDegree, final int localIndex)
	{
		shapeFunctionsOnCell = new HashMap<>();
		faces = new HashSet<>();
		this.polynomialDegree = polynomialDegree;
		firstDefinedFunction = TPShapeFunction.ReferenceFunction(cell.getDimension(), polynomialDegree,
		                                                         localIndex);
		new TPShapeFunction(cell.getReferenceCell(), polynomialDegree, localIndex);
		firstDefinedCell = cell;
		functional =
			new LagrangeNodeFunctional(firstDefinedCell.transformFromReferenceCell(firstDefinedFunction.nodeFunctional.getPoint()));
		shapeFunctionsOnCell.put(cell, firstDefinedFunction);
		addFaces(nodeOnWhichFaces(cell));
		shapeFunctionsOnCell.keySet()
		                    .forEach(c -> faces.addAll(c.getFaces()));
		cellCodes = new int[shapeFunctionsOnCell.size()];
		shapeFunctionsOnCellFast = new TPShapeFunction[cellCodes.length];
		int i = 0;
		for (final var sf : shapeFunctionsOnCell.entrySet())
		{
			cellCodes[i] = sf.getKey()
			                 .doneCode();
			shapeFunctionsOnCellFast[i] = sf.getValue();
			i++;
		}
	}
	
	private List<TPFace> nodeOnWhichFaces(final TPCell cell)
	{
		return cell.getFaces()
		           .stream()
		           .filter(getNodeFunctional()::usesFace)
		           .collect(Collectors.toList());
	}
	
	private void addFaces(final List<TPFace> newFaces)
	{
		for (final TPFace face : newFaces)
		{
			if (faces.contains(face))
				continue;
			faces.add(face);
			for (final TPCell cellOfFace : face.getCells())
			{
				if (shapeFunctionsOnCell.containsKey(cellOfFace))
					continue;
				final TPShapeFunction newFunction =
					TPShapeFunction
						.ReferenceFunction(
							cellOfFace.getDimension(),
							polynomialDegree,
							TPShapeFunction.getIndexFromReferenceNodePoint(
								polynomialDegree,
								cellOfFace.transformToReferenceCell(
									getNodeFunctional().getPoint())));
				shapeFunctionsOnCell.put(cellOfFace, newFunction);
				addFaces(nodeOnWhichFaces(cellOfFace));
			}
		}
	}
	
	TPShapeFunction getFunction(final TPCell cell)
	{
		final int cellcode = cell.doneCode();
		for (int i = 0; i < cellCodes.length; i++)
		{
			if (cellcode == cellCodes[i])
				return shapeFunctionsOnCellFast[i];
		}
		return null;
	}
	
	@Override
	public Double valueOnReferenceCell(final CoordinateVector pos, final TPCell cell)
	{
		if (PerformanceArguments.getInstance().executeChecks)
			if (!cell.getReferenceCell()
			         .isInCell(pos)) throw new IllegalArgumentException("pos is not in cell");
		return getFunction(cell)
			.fastValueOnReferenceCell(pos);
	}
	
	public LagrangeNodeFunctional nodeFunctionalOnReferenceCell(final TPCell cell)
	{
		return getFunction(cell)
			.getNodeFunctional();
	}
	
	@Override
	public CoordinateVector gradientOnReferenceCell(final CoordinateVector pos, final TPCell cell)
	{
		final double[] vals = getFunction(cell)
			.fastGradientOnReferenceCell(pos);
		for (int i = 0; i < vals.length; i++)
			vals[i] *= 1. / cell.size[i];
		return new CoordinateVector(vals, true);
	}
	
	@Override
	public List<TPCell> getCells()
	{
		return new ArrayList<>(shapeFunctionsOnCell.keySet());
	}
	
	@Override
	public Set<TPFace> getFaces()
	{
		return faces;
	}
	
	@Override
	public LagrangeNodeFunctional getNodeFunctional()
	{
		return functional;
	}
	
	@Override
	public int getGlobalIndex()
	{
		return globalIndex;
	}
	
	public void setGlobalIndex(final int index)
	{
		globalIndex = index;
	}
	
	@Override
	public boolean equals(final Object o)
	{
		if (this == o) return true;
		if (o == null || getClass() != o.getClass()) return false;
		final CTPShapeFunction that = (CTPShapeFunction) o;
		return Objects.equals(getNodeFunctional().getPoint(),
		                      that.getNodeFunctional()
		                          .getPoint());
	}
	
	@Override
	public int compareTo(@NotNull final CTPShapeFunction o)
	{
		return getNodeFunctional().getPoint()
		                          .compareTo(o.getNodeFunctional()
		                                      .getPoint());
	}
	
	@Override
	public int hashCode()
	{
		return getNodeFunctional().getPoint()
		                          .hashCode();
	}
	
	@Override
	public String toString()
	{
		return "CTPShapeFunction{ " + "functionalpoint:  " + getNodeFunctional()
			.getPoint() + '}';
	}
}
