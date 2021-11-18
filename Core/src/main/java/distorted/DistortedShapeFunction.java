package distorted;

import basic.DoubleCompare;
import basic.LagrangeNodeFunctional;
import basic.PerformanceArguments;
import basic.ScalarShapeFunction;
import distorted.geometry.DistortedCell;
import distorted.geometry.DistortedFace;
import it.unimi.dsi.fastutil.ints.Int2ObjectArrayMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectMap;
import linalg.CoordinateMatrix;
import linalg.CoordinateVector;
import org.jetbrains.annotations.NotNull;
import tensorproduct.TPShapeFunction;
import tensorproduct.geometry.TPFace;

import java.util.*;

public class DistortedShapeFunction
	implements ScalarShapeFunction<DistortedCell, DistortedFace>,
	DistortedScalarFunction,
	Comparable<DistortedShapeFunction>

{
	private final Map<DistortedCell, TPShapeFunction> shapeFunctionsOnCell;
	private final Int2ObjectMap<TPShapeFunction> shapeFunctionsOnCellFast;
	private final Set<DistortedFace> faces;
	private final TPShapeFunction firstDefinedFunction;
	private final DistortedCell firstDefinedCell;
	int polynomialDegree;
	private int globalIndex;
	
	public DistortedShapeFunction(final DistortedCell cell, final int polynomialDegree, final int localIndex)
	{
		shapeFunctionsOnCell = new HashMap<>();
		shapeFunctionsOnCellFast = new Int2ObjectArrayMap<>();
		faces = new HashSet<>();
		this.polynomialDegree = polynomialDegree;
		firstDefinedFunction = new TPShapeFunction(cell.referenceCell, polynomialDegree, localIndex);
		firstDefinedCell = cell;
		shapeFunctionsOnCell.put(cell, firstDefinedFunction);
		shapeFunctionsOnCellFast.put(cell.doneCode(), firstDefinedFunction);
		checkIfPointOnFace(cell);
	}
	
	private void checkIfPointOnFace(final DistortedCell cell)
	{
		
		for (final DistortedFace face : cell.getFaces())
		{
			final TPFace referenceFace = TPFace.fromVertices(cell.getReferenceVerticesOfFace(face),
			                                                 face.isBoundaryFace());
			if (shapeFunctionsOnCell.get(cell)
			                        .getNodeFunctional()
			                        .usesFace(referenceFace))
			{
				if (faces.add(face))
				{
					for (final DistortedCell cellOfFace : face.getCells())
					{
						final TPShapeFunction newFunction = new TPShapeFunction(
							cell.referenceCell, polynomialDegree,
							cellOfFace.transformPreciseToReferenceCell(
								getNodeFunctional().getPoint()));
						shapeFunctionsOnCell.put(cellOfFace, newFunction);
						shapeFunctionsOnCellFast.put(cellOfFace.doneCode(), newFunction);
						checkIfPointOnFace(cellOfFace);
					}
				}
			}
		}
	}
	
	@Override
	public Double valueOnReferenceCell(final CoordinateVector pos, final DistortedCell cell)
	{
		if (PerformanceArguments.getInstance().executeChecks)
			if (!cell.referenceCell.isInCell(pos)) throw new IllegalArgumentException("pos is not in cell");
		//if (shapeFunctionsOnCell.containsKey(cell))
		return shapeFunctionsOnCellFast.get(cell.doneCode())
		                               .fastValueInCell(pos);
		//return 0;
	}
	
	public LagrangeNodeFunctional nodeFunctionalOnReferenceCell(final DistortedCell cell)
	{
		return shapeFunctionsOnCellFast.get(cell.doneCode())
		                               .getNodeFunctional();
	}
	
	@Override
	public CoordinateVector gradientOnReferenceCell(final CoordinateVector pos, final DistortedCell cell)
	{
		if (PerformanceArguments.getInstance().executeChecks)
			if (!cell.referenceCell.isInCell(pos)) throw new IllegalArgumentException("pos is not in cell");
		//if (shapeFunctionsOnCell.containsKey(cell))
		//{
		final CoordinateMatrix mat = cell.transformationGradientFromReferenceCell(pos)
		                                 .inverse();
		final CoordinateVector grad = shapeFunctionsOnCellFast.get(cell.doneCode())
		                                                      .gradient(pos);
		return mat.mvMul(grad);
		//}
		//return new CoordinateVector(pos.getLength());
	}
	
	@Override
	public Set<DistortedCell> getCells()
	{
		return shapeFunctionsOnCell.keySet();
	}
	
	@Override
	public Set<DistortedFace> getFaces()
	{
		return faces;
	}
	
	@Override
	public LagrangeNodeFunctional getNodeFunctional()
	{
		return new LagrangeNodeFunctional(firstDefinedCell.transformFromReferenceCell(
			firstDefinedFunction.getNodeFunctional()
			                    .getPoint()));
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
		final DistortedShapeFunction that = (DistortedShapeFunction) o;
		return Objects.equals(getNodeFunctional().getPoint(),
		                      that.getNodeFunctional()
		                          .getPoint());
	}
	
	@Override
	public int compareTo(@NotNull final DistortedShapeFunction o)
	{
		return getNodeFunctional().getPoint()
		                          .compareTo(o.getNodeFunctional()
		                                      .getPoint());
	}
	
	@Override
	public String toString()
	{
		return "DistortedShapeFunction{" + "functionalpoint: " + getNodeFunctional()
			.getPoint() + " ::: \t" + getNodeFunctional().getPoint()
		                                                     .at(0) + " " + DoubleCompare.doubleHash(
			getNodeFunctional().getPoint()
			                   .at(0)) + " \t" + getNodeFunctional()
			.getPoint()
			.at(1) + " " + DoubleCompare.doubleHash(
			getNodeFunctional().getPoint()
			                   .at(1)) + " \t" + hashCode() + '}' + "\n";
	}
}
