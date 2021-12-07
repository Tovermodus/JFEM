package distorted;

import basic.AcceptsMatrixBoundaryValues;
import basic.Assembleable;
import basic.ShapeFunction;
import distorted.geometry.CircleGrid;
import distorted.geometry.DistortedCell;
import distorted.geometry.DistortedFace;
import linalg.CoordinateVector;

public abstract class CircleGridSpace<ST extends ShapeFunction<DistortedCell, DistortedFace, valueT, gradientT, hessianT> & Comparable<ST>, valueT, gradientT, hessianT>
	extends DistortedGridSpace<CircleGrid, ST, valueT, gradientT, hessianT>
	implements AcceptsMatrixBoundaryValues<DistortedCell, DistortedFace, ST, valueT, gradientT, hessianT>,
	Assembleable

{
	public final double radius;
	public final CoordinateVector center;
	
	public CircleGridSpace(final CoordinateVector center, final double radius, final int refinements)
	{
		super(new CircleGrid(center, radius, refinements));
		this.radius = radius;
		this.center = center;
	}
}
