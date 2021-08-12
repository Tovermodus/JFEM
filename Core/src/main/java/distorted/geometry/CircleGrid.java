package distorted.geometry;

import basic.DoubleCompare;
import com.google.common.collect.ImmutableSet;
import linalg.CoordinateVector;
import linalg.IntCoordinates;

import java.util.*;
import java.util.stream.Collectors;

public class CircleGrid
{
	public final int dimension;
	private final CoordinateVector centerPoint;
	private final double radius;
	public ImmutableSet<DistortedCell> cells;
	public ImmutableSet<DistortedFace> faces;
	
	public CircleGrid(final int dimension, final CoordinateVector centerPoint, final double radius, final int refinements)
	{
		this.dimension = dimension;
		this.centerPoint = centerPoint;
		this.radius = radius;
		final HashSet<DistortedCell> genCells = generateCells();
		final HashSet<DistortedFace> genFaces = generateFaces(genCells);
		for (int i = 0; i < refinements; i++)
			refine(genCells, genFaces);
		cells = ImmutableSet.copyOf(genCells);
		faces = ImmutableSet.copyOf(genFaces);
	}
	
	public HashSet<DistortedCell> generateCells()
	{
		if (dimension == 2)
		{
			final double sideLength = radius * 0.5 / Math.sqrt(2);
			final CoordinateVector center_ll = CoordinateVector.fromValues(-sideLength, -sideLength);
			final CoordinateVector center_lr = CoordinateVector.fromValues(sideLength, -sideLength);
			final CoordinateVector center_ur = CoordinateVector.fromValues(sideLength, sideLength);
			final CoordinateVector center_ul = CoordinateVector.fromValues(-sideLength, sideLength);
			final CoordinateVector outer_ll = center_ll.mul(2);
			final CoordinateVector outer_lr = center_lr.mul(2);
			final CoordinateVector outer_ur = center_ur.mul(2);
			final CoordinateVector outer_ul = center_ul.mul(2);
			center_ll.addInPlace(centerPoint);
			center_lr.addInPlace(centerPoint);
			center_ur.addInPlace(centerPoint);
			center_ul.addInPlace(centerPoint);
			outer_ll.addInPlace(centerPoint);
			outer_lr.addInPlace(centerPoint);
			outer_ur.addInPlace(centerPoint);
			outer_ul.addInPlace(centerPoint);
			final DistortedCell center =
				new DistortedCell(List
					                  .of(center_ll, center_lr, center_ur, center_ul)
					                  .toArray(new CoordinateVector[4]));
			final DistortedCell left =
				new DistortedCell(List
					                  .of(outer_ll, center_ll, center_ul, outer_ul)
					                  .toArray(new CoordinateVector[4]));
			final DistortedCell right =
				new DistortedCell(List
					                  .of(center_lr, outer_lr, outer_ur, center_ur)
					                  .toArray(new CoordinateVector[4]));
			final DistortedCell top =
				new DistortedCell(List
					                  .of(center_ul, center_ur, outer_ur, outer_ul)
					                  .toArray(new CoordinateVector[4]));
			final DistortedCell bottom =
				new DistortedCell(List
					                  .of(outer_ll, outer_lr, center_lr, center_ll)
					                  .toArray(new CoordinateVector[4]));
			return new HashSet<>(List.of(center, left, right, top, bottom));
		} else if (dimension == 3)
		{
			final double sideLength = radius * 0.5 / Math.sqrt(3);
			final CoordinateVector center_llf = CoordinateVector.fromValues(-sideLength, -sideLength,
			                                                                -sideLength);
			final CoordinateVector center_lrf = CoordinateVector.fromValues(sideLength, -sideLength, -sideLength);
			final CoordinateVector center_lrb = CoordinateVector.fromValues(sideLength, sideLength, -sideLength);
			final CoordinateVector center_llb = CoordinateVector.fromValues(-sideLength, sideLength, -sideLength);
			final CoordinateVector center_ulf = CoordinateVector.fromValues(-sideLength, -sideLength, sideLength);
			final CoordinateVector center_urf = CoordinateVector.fromValues(sideLength, -sideLength, sideLength);
			final CoordinateVector center_urb = CoordinateVector.fromValues(sideLength, sideLength, sideLength);
			final CoordinateVector center_ulb = CoordinateVector.fromValues(-sideLength, sideLength, sideLength);
			final CoordinateVector outer_llf = center_llf.mul(2);
			final CoordinateVector outer_lrf = center_lrf.mul(2);
			final CoordinateVector outer_urf = center_urf.mul(2);
			final CoordinateVector outer_ulf = center_ulf.mul(2);
			final CoordinateVector outer_llb = center_llb.mul(2);
			final CoordinateVector outer_lrb = center_lrb.mul(2);
			final CoordinateVector outer_urb = center_urb.mul(2);
			final CoordinateVector outer_ulb = center_ulb.mul(2);
			center_llf.addInPlace(centerPoint);
			center_lrf.addInPlace(centerPoint);
			center_lrb.addInPlace(centerPoint);
			center_llb.addInPlace(centerPoint);
			center_ulf.addInPlace(centerPoint);
			center_urf.addInPlace(centerPoint);
			center_urb.addInPlace(centerPoint);
			center_ulb.addInPlace(centerPoint);
			outer_llf.addInPlace(centerPoint);
			outer_lrf.addInPlace(centerPoint);
			outer_urf.addInPlace(centerPoint);
			outer_ulf.addInPlace(centerPoint);
			outer_llb.addInPlace(centerPoint);
			outer_lrb.addInPlace(centerPoint);
			outer_urb.addInPlace(centerPoint);
			outer_ulb.addInPlace(centerPoint);
			final DistortedCell center =
				new DistortedCell(List.of(center_llf, center_lrf, center_lrb, center_llb,
				                          center_ulf, center_urf, center_urb, center_ulb)
				                      .toArray(new CoordinateVector[8]));
			final DistortedCell left =
				new DistortedCell(List.of(outer_llf, center_llf, center_llb, outer_llb,
				                          outer_ulf, center_ulf, center_ulb, outer_ulb)
				                      .toArray(new CoordinateVector[8]));
			final DistortedCell back =
				new DistortedCell(List.of(outer_llb, center_llb, center_lrb, outer_lrb
					, outer_ulb, center_ulb, center_urb, outer_urb)
				                      .toArray(new CoordinateVector[8]));
			final DistortedCell bottom =
				new DistortedCell(List.of(outer_llf, outer_lrf, outer_lrb, outer_llb,
				                          center_llf, center_lrf, center_lrb, center_llb)
				                      .toArray(new CoordinateVector[8]));
			final DistortedCell front =
				new DistortedCell(List.of(outer_llf, outer_lrf, center_lrf, center_llf,
				                          outer_ulf, outer_urf, center_urf, center_ulf)
				                      .toArray(new CoordinateVector[8]));
			final DistortedCell right =
				new DistortedCell(List.of(outer_lrf, outer_lrb, center_lrb, center_lrf,
				                          outer_urf, outer_urb, center_urb, center_urf)
				                      .toArray(new CoordinateVector[8]));
			final DistortedCell top =
				new DistortedCell(List.of(center_ulf, center_urf, center_urb, center_ulb,
				                          outer_ulf, outer_urf, outer_urb, outer_ulb)
				                      .toArray(new CoordinateVector[8]));
			return new HashSet<>(List.of(center, bottom, top, front, right, back, left));
		} else throw new IllegalArgumentException("Only supports dimensions 2 and 3");
	}
	
	public HashSet<DistortedFace> generateFaces(final HashSet<DistortedCell> genCells)
	{
		final HashSet<DistortedFace> faces = new HashSet<>();
		for (final DistortedCell cell : genCells)
		{
			for (final CoordinateVector[] vertices : cell.getSides())
			{
				final boolean isOnBoundary = Arrays.stream(vertices).mapToDouble(
					CoordinateVector::euclidianNorm)
				                                   .allMatch(norm -> DoubleCompare.almostEqual(norm, radius));
				final DistortedFace face = new DistortedFace(vertices, isOnBoundary);
				for (final DistortedCell c : getCellsBelongingToFace(face, genCells))
				{
					cell.faces.add(face);
					face.addCell(c);
				}
				faces.add(face);
			}
		}
		return faces;
	}
	
	private HashSet<DistortedCell> getCellsBelongingToFace(final DistortedFace f, final Collection<DistortedCell> from)
	{
		final HashSet<DistortedCell> cellsBelonging = new HashSet<>();
		for (final DistortedCell cell : from)
		{
			if (f.isOnCell(cell))
				cellsBelonging.add(cell);
		}
		return cellsBelonging;
	}
	
	public void refine(final HashSet<DistortedCell> genCells, final HashSet<DistortedFace> genFaces)
	{
//		wennamrand erstrandkoordinate bestimmen
//		für mitte alles halbieren
//
//		für faces: naiv wie oben aber nur faces von den cells die vorher über die face verbunden waren
	}
	
	public List<CoordinateVector> generatePlotPoints(final int resolution)
	{
		final int totalPoints = (int) Math.pow(resolution, dimension);
		final int pointsPerCell = totalPoints / cells.size();
		final int pointsPerDirectionPerCell = (int) (Math.pow(pointsPerCell, 1. / dimension));
		final List<CoordinateVector> ret = new ArrayList<>();
		final IntCoordinates pointsPerDimension = IntCoordinates.repeat(pointsPerDirectionPerCell, dimension);
		for (final DistortedCell cell : cells)
		{
			ret.addAll(pointsPerDimension
				           .range()
				           .stream()
				           .map(c ->
				                {
					                final CoordinateVector coords = new CoordinateVector(dimension);
					                for (int i = 0; i < dimension; i++)
					                {
						                coords.set(1. / (pointsPerDirectionPerCell - 1) *
							                           (c.get(i)), i);
					                }
					                return coords;
				                })
				           .map(cell::transformFromReferenceCell)
				           .collect(Collectors.toList()));
		}
		return ret;
	}
}
