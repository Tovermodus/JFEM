package distorted;

import basic.DoubleCompare;
import basic.FastEvaluatedScalarShapeFunction;
import basic.LagrangeNodeFunctional;
import basic.PerformanceArguments;
import distorted.geometry.DistortedCell;
import distorted.geometry.DistortedFace;
import linalg.CoordinateMatrix;
import linalg.CoordinateVector;
import org.jetbrains.annotations.NotNull;
import tensorproduct.TPShapeFunction;
import tensorproduct.geometry.TPFace;

import java.util.*;

public class DistortedShapeFunction implements FastEvaluatedScalarShapeFunction<DistortedCell, DistortedFace>, Comparable<DistortedShapeFunction>

{
	private final Map<DistortedCell, TPShapeFunction> shapeFunctionsOnCell;
	private final Set<DistortedFace> faces;
	private final TPShapeFunction firstDefinedFunction;
	private final DistortedCell firstDefinedCell;
	int polynomialDegree;
	private int globalIndex;
	
	public DistortedShapeFunction(final DistortedCell cell, final int polynomialDegree, final int localIndex)
	{
		shapeFunctionsOnCell = new HashMap<>();
		faces = new HashSet<>();
		this.polynomialDegree = polynomialDegree;
		firstDefinedFunction = new TPShapeFunction(cell.referenceCell, polynomialDegree, localIndex);
		firstDefinedCell = cell;
		shapeFunctionsOnCell.put(cell, firstDefinedFunction);
		checkIfPointOnFace(cell);
	}
	
	private void checkIfPointOnFace(final DistortedCell cell)
	{
		
		for (final DistortedFace face : cell.getFaces())
		{
			final TPFace referenceFace = TPFace.fromVertices(cell.getReferenceVerticesOfFace(face),
			                                                 face.isBoundaryFace());
			if (shapeFunctionsOnCell.get(cell).getNodeFunctional().usesFace(referenceFace))
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
						checkIfPointOnFace(cellOfFace);
					}
				}
			}
		}
	}
	
	@Override
	public double fastValueInCell(final CoordinateVector pos, final DistortedCell cell)
	{
		return valueInCell(pos, cell);
	}
	
	@Override
	public double[] fastGradientInCell(final CoordinateVector pos, final DistortedCell cell)
	{
		return gradientInCell(pos, cell).toArray();
	}
	
	@Override
	public Double valueInCell(final CoordinateVector pos, final DistortedCell cell)
	{
		return valueOnReferenceCell(cell.transformToReferenceCell(pos), cell);
	}
	
	public double valueOnReferenceCell(final CoordinateVector pos, final DistortedCell cell)
	{
		if (PerformanceArguments.getInstance().executeChecks)
			if (!cell.referenceCell.isInCell(pos)) throw new IllegalArgumentException("pos is not in cell");
		//if (shapeFunctionsOnCell.containsKey(cell))
		return shapeFunctionsOnCell.get(cell).fastValue(pos);
		//return 0;
	}
	
	public LagrangeNodeFunctional nodeFunctionalOnReferenceCell(final DistortedCell cell)
	{
		return shapeFunctionsOnCell.get(cell).getNodeFunctional();
	}
	
	public CoordinateVector gradientOnReferenceCell(final CoordinateVector pos, final DistortedCell cell)
	{
		if (PerformanceArguments.getInstance().executeChecks)
			if (!cell.referenceCell.isInCell(pos)) throw new IllegalArgumentException("pos is not in cell");
		//if (shapeFunctionsOnCell.containsKey(cell))
		//{
		final CoordinateMatrix mat = cell.transformationGradientFromReferenceCell(pos).inverse();
		final CoordinateVector grad = shapeFunctionsOnCell.get(cell).gradient(pos);
		return mat.tvMul(grad);
		//}
		//return new CoordinateVector(pos.getLength());
	}
	
	@Override
	public CoordinateVector gradientInCell(final CoordinateVector pos, final DistortedCell cell)
	{
		return gradientOnReferenceCell(cell.transformToReferenceCell(pos), cell);
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
			firstDefinedFunction.getNodeFunctional().getPoint()));
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
		return Objects.equals(getNodeFunctional().getPoint(), that.getNodeFunctional().getPoint());
	}
	
	@Override
	public int compareTo(@NotNull final DistortedShapeFunction o)
	{
		return getNodeFunctional().getPoint().compareTo(o.getNodeFunctional().getPoint());
	}
	
	@Override
	public String toString()
	{
		return "DistortedShapeFunction{" + "functionalpoint: " + getNodeFunctional()
			.getPoint() + " ::: \t" + getNodeFunctional().getPoint().at(0) + " " + DoubleCompare.doubleHash(
			getNodeFunctional().getPoint().at(0)) + " \t" + getNodeFunctional()
			.getPoint()
			.at(1) + " " + DoubleCompare.doubleHash(
			getNodeFunctional().getPoint().at(1)) + " \t" + hashCode() + '}' + "\n";
	}
}
