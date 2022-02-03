package distorted;

import basic.AcceptsMatrixBoundaryValues;
import basic.Assembleable;
import basic.ShapeFunction;
import distorted.geometry.BrickGrid;
import distorted.geometry.DistortedCell;
import distorted.geometry.DistortedFace;
import linalg.CoordinateVector;

public abstract class BrickGridSpace<ST extends ShapeFunction<DistortedCell, DistortedFace, valueT, gradientT, hessianT> & Comparable<ST>, valueT, gradientT, hessianT>
	extends DistortedGridSpace<ST, valueT, gradientT, hessianT>
	implements AcceptsMatrixBoundaryValues<DistortedCell, DistortedFace, ST, valueT, gradientT, hessianT>,
	Assembleable

{
	public final double width;
	public final double height;
	public final CoordinateVector center;
	
	public BrickGridSpace(final CoordinateVector center, final double width, final double height,
	                      final int refinements)
	{
		super(new BrickGrid(center, width, height, refinements));
		this.width = width;
		this.height = height;
		this.center = center;
	}
}
