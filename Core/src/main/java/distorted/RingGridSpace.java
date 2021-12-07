package distorted;

import basic.AcceptsMatrixBoundaryValues;
import basic.Assembleable;
import basic.ShapeFunction;
import distorted.geometry.DistortedCell;
import distorted.geometry.DistortedFace;
import distorted.geometry.RingGrid;
import linalg.CoordinateVector;

public abstract class RingGridSpace<ST extends ShapeFunction<DistortedCell, DistortedFace, valueT, gradientT, hessianT> & Comparable<ST>, valueT, gradientT, hessianT>
	extends DistortedGridSpace<RingGrid, ST, valueT, gradientT, hessianT>
	implements AcceptsMatrixBoundaryValues<DistortedCell, DistortedFace, ST, valueT, gradientT, hessianT>,
	Assembleable

{
	public final double innerRadius;
	public final double outerRadius;
	public final CoordinateVector center;
	
	public RingGridSpace(final CoordinateVector center, final double innerRadius, final double outerRadius,
	                     final int refinements)
	{
		super(new RingGrid(center, innerRadius, outerRadius, refinements));
		this.innerRadius = innerRadius;
		this.outerRadius = outerRadius;
		this.center = center;
	}
}
