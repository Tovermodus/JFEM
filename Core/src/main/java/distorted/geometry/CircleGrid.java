package distorted.geometry;

import com.google.common.collect.ImmutableList;
import io.vavr.collection.Array;
import linalg.CoordinateVector;
import linalg.IntCoordinates;
import linalg.Matrix;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

public class CircleGrid
{
	public ImmutableList<DistortedCell> cells;
	public ImmutableList<DistortedFace> faces;
	public final int dimension;
	private final CoordinateVector center;
	private final double radius;
	
	public CircleGrid(int dimension, CoordinateVector center, double radius, int refinements)
	{
		this.dimension = dimension;
		this.center = center;
		this.radius = radius;
		List<DistortedCell> genCells = generateCells();
		List<DistortedFace> genFaces = generateFaces(genCells);
		for(int i = 0; i < refinements; i++)
			refine(genCells, genFaces);
		cells = ImmutableList.copyOf(genCells);
		faces = ImmutableList.copyOf(genFaces);
	}
	public List<DistortedCell> generateCells()
	{
		if(dimension == 2)
		{
			double sideLength = radius*0.5/Math.sqrt(2);
			CoordinateVector center_ll = CoordinateVector.fromValues(-sideLength,-sideLength);
			CoordinateVector center_lr = CoordinateVector.fromValues(sideLength,-sideLength);
			CoordinateVector center_ur = CoordinateVector.fromValues(sideLength,sideLength);
			CoordinateVector center_ul = CoordinateVector.fromValues(-sideLength,sideLength);
			CoordinateVector outer_ll = center_ll.mul(2);
			CoordinateVector outer_lr = center_lr.mul(2);
			CoordinateVector outer_ur = center_ur.mul(2);
			CoordinateVector outer_ul = center_ul.mul(2);
			DistortedCell center =
				new DistortedCell(List.of(center_ll, center_lr, center_ur, center_ul).toArray(new CoordinateVector[4]));
			DistortedCell left =
				new DistortedCell(List.of(outer_ll, center_ll, center_ul, outer_ul).toArray(new CoordinateVector[4]));
			DistortedCell right =
				new DistortedCell(List.of(center_lr, outer_lr, outer_ur, center_ur).toArray(new CoordinateVector[4]));
			DistortedCell top =
				new DistortedCell(List.of(center_ul, center_ur, outer_ur, outer_ul).toArray(new CoordinateVector[4]));
			DistortedCell bottom =
				new DistortedCell(List.of(outer_ll, outer_lr, center_lr, center_ll).toArray(new CoordinateVector[4]));
			return List.of(center, left, right, top, bottom);
		}
		else if(dimension == 3)
		{
			double sideLength = radius*0.5/Math.sqrt(3);
			CoordinateVector center_llf = CoordinateVector.fromValues(-sideLength,-sideLength,-sideLength);
			CoordinateVector center_lrf = CoordinateVector.fromValues(sideLength,-sideLength,-sideLength);
			CoordinateVector center_lrb = CoordinateVector.fromValues(sideLength,sideLength,-sideLength);
			CoordinateVector center_llb = CoordinateVector.fromValues(-sideLength,sideLength,-sideLength);
			CoordinateVector center_ulf = CoordinateVector.fromValues(-sideLength,-sideLength,sideLength);
			CoordinateVector center_urf = CoordinateVector.fromValues(sideLength,-sideLength,sideLength);
			CoordinateVector center_urb = CoordinateVector.fromValues(sideLength,sideLength,sideLength);
			CoordinateVector center_ulb = CoordinateVector.fromValues(-sideLength,sideLength,sideLength);
			CoordinateVector outer_llf = center_llf.mul(2);
			CoordinateVector outer_lrf = center_lrf.mul(2);
			CoordinateVector outer_urf = center_urf.mul(2);
			CoordinateVector outer_ulf = center_ulf.mul(2);
			CoordinateVector outer_llb = center_llb.mul(2);
			CoordinateVector outer_lrb = center_lrb.mul(2);
			CoordinateVector outer_urb = center_urb.mul(2);
			CoordinateVector outer_ulb = center_ulb.mul(2);
			DistortedCell center =
				new DistortedCell(List.of(center_llf, center_lrf, center_lrb, center_llb,
							center_ulf,center_urf, center_urb, center_ulb)
					.toArray(new CoordinateVector[8]));
			DistortedCell left =
				new DistortedCell(List.of(outer_llf, center_llf, center_llb, outer_llb,
					outer_ulf, center_ulf, center_ulb, outer_ulb)
					.toArray(new CoordinateVector[8]));
			DistortedCell back =
				new DistortedCell(List.of(  outer_llb,center_llb, center_lrb, outer_lrb
					, outer_ulb, center_ulb, center_urb, outer_urb)
					.toArray(new CoordinateVector[8]));
			DistortedCell bottom =
				new DistortedCell(List.of(outer_llf, outer_lrf, outer_lrb, outer_llb,
					center_llf,center_lrf, center_lrb, center_llb)
					.toArray(new CoordinateVector[8]));
			DistortedCell front =
				new DistortedCell(List.of(outer_llf, outer_lrf, center_lrf, center_llf,
					outer_ulf,outer_urf, center_urf, center_ulf)
					.toArray(new CoordinateVector[8]));
			DistortedCell right =
				new DistortedCell(List.of(outer_lrf, outer_lrb, center_lrb, center_lrf,
					outer_urf, outer_urb, center_urb,  center_urf)
					.toArray(new CoordinateVector[8]));
			DistortedCell top =
				new DistortedCell(List.of(center_ulf, center_urf, center_urb, center_ulb,
					outer_ulf, outer_urf, outer_urb, outer_ulb)
					.toArray(new CoordinateVector[8]));
			return List.of(center, bottom, top, front, right,back, left);
		}
		else throw new IllegalArgumentException("Only supports dimensions 2 and 3");
	}
	public List<DistortedFace> generateFaces( List<DistortedCell> genCells)
	{
		return new ArrayList<>();
	}
	public void refine( List<DistortedCell> genCells, List<DistortedFace> genFaces)
	{
	
	}
	
	public List<CoordinateVector> generatePlotPoints(int resolution)
	{
		int totalPoints = (int)Math.pow(resolution, dimension);
		int pointsPerCell = totalPoints;
		int pointsPerDirectionPerCell = (int)(Math.pow(pointsPerCell, 1./dimension));
		List<CoordinateVector> ret = new ArrayList<>();
		IntCoordinates pointsPerDimension = IntCoordinates.repeat(pointsPerDirectionPerCell, dimension);
		for(DistortedCell cell: cells)
		{
			ret.addAll(pointsPerDimension.range().stream().map(c ->
			{
				CoordinateVector coords = new CoordinateVector(dimension);
				for (int i = 0; i < dimension; i++)
				{
					coords.set(1. / (pointsPerDirectionPerCell - 1) * (c.get(i)), i);
				}
				return coords;
			}).map(cell::transformFromReferenceCell).collect(Collectors.toList()));
		}
		return ret;
	}
}
