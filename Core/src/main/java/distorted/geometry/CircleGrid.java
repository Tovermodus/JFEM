package distorted.geometry;

import basic.DoubleCompare;
import com.google.common.collect.HashMultimap;
import com.google.common.collect.ImmutableSet;
import com.google.common.collect.Multimap;
import linalg.CoordinateVector;
import linalg.IntCoordinates;
import org.jetbrains.annotations.NotNull;

import java.util.*;
import java.util.stream.Collectors;

public class CircleGrid
{
	public final int dimension;
	private final CoordinateVector centerPoint;
	private final double radius;
	public ImmutableSet<DistortedCell> cells;
	public ImmutableSet<DistortedFace> faces;
	
	public CircleGrid(final CoordinateVector centerPoint, final double radius, final int refinements)
	{
		this.dimension = centerPoint.getLength();
		this.centerPoint = centerPoint;
		this.radius = radius;
		
		Multimap<DistortedCell, DistortedCell> refinedCells = generateCells();
		Collection<DistortedFace> genFaces = generateFaces(refinedCells);
		for (int i = 0; i < refinements; i++)
		{
			System.out.println("refine");
			refinedCells = refineCells(refinedCells);
			genFaces = refineFaces(refinedCells);
			for (final DistortedFace f : genFaces)
			{
				if (!f.isBoundaryFace() && f.getCells().size() != 2)
					System.out.println(f);
			}
		}
		cells = ImmutableSet.copyOf(refinedCells.values());
		faces = ImmutableSet.copyOf(genFaces);
	}
	
	private static HashSet<DistortedCell> getCellsBelongingToFace(final DistortedFace f, final Collection<DistortedCell> from)
	{
		final HashSet<DistortedCell> cellsBelonging = new HashSet<>();
		for (final DistortedCell cell : from)
		{
			if (f.isOnCell(cell))
				cellsBelonging.add(cell);
		}
		return cellsBelonging;
	}
	
	public Multimap<DistortedCell, DistortedCell> generateCells()
	{
		final Multimap<DistortedCell, DistortedCell> ret = HashMultimap.create();
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
			ret.put(center, center);
			ret.put(left, left);
			ret.put(right, right);
			ret.put(top, top);
			ret.put(bottom, bottom);
			return ret;
		} else if (dimension == 3)
		{
			final double sideLength = radius * 0.5 / Math.sqrt(3);
			final CoordinateVector center_llf = CoordinateVector.fromValues(-sideLength, -sideLength,
			                                                                -sideLength);
			final CoordinateVector center_lrf = CoordinateVector.fromValues(sideLength, -sideLength,
			                                                                -sideLength);
			final CoordinateVector center_lrb = CoordinateVector.fromValues(sideLength, sideLength,
			                                                                -sideLength);
			final CoordinateVector center_llb = CoordinateVector.fromValues(-sideLength, sideLength,
			                                                                -sideLength);
			final CoordinateVector center_ulf = CoordinateVector.fromValues(-sideLength, -sideLength,
			                                                                sideLength);
			final CoordinateVector center_urf = CoordinateVector.fromValues(sideLength, -sideLength,
			                                                                sideLength);
			final CoordinateVector center_urb = CoordinateVector.fromValues(sideLength, sideLength,
			                                                                sideLength);
			final CoordinateVector center_ulb = CoordinateVector.fromValues(-sideLength, sideLength,
			                                                                sideLength);
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
			ret.put(center, center);
			ret.put(left, left);
			ret.put(back, back);
			ret.put(bottom, bottom);
			ret.put(front, front);
			ret.put(right, right);
			ret.put(top, top);
			return ret;
		} else throw new IllegalArgumentException("Only supports dimensions 2 and 3");
	}
	
	public Collection<DistortedFace> generateFaces(final Multimap<DistortedCell, DistortedCell> genCells)
	{
		final HashSet<DistortedFace> faces = new HashSet<>();
		for (final DistortedCell cell : genCells.values())
		{
			createFaces(faces, cell, genCells.values());
		}
		return faces;
	}
	
	public HashSet<DistortedFace> refineFaces(final Multimap<DistortedCell, DistortedCell> refinedCells)
	{
		final HashSet<DistortedFace> faces = new HashSet<>();
		for (final DistortedCell cell : refinedCells.keySet())
		{
			final List<DistortedCell> refinedNeighbouringCells = getRefinedNeighbouringCells(refinedCells,
			                                                                                 cell);
			for (final DistortedCell refinedCell : refinedCells.get(cell))
			{
				createFaces(faces, refinedCell, refinedNeighbouringCells);
			}
		}
		return faces;
	}
	
	private void createFaces(final HashSet<DistortedFace> faces, final DistortedCell cell,
	                         final Collection<DistortedCell> neighbouringCells)
	{
		createFaces(faces, cell, neighbouringCells, radius, centerPoint);
	}
	
	public static void createFaces(final HashSet<DistortedFace> faces, final DistortedCell cell,
	                               final Collection<DistortedCell> neighbouringCells, final double radius,
	                               final CoordinateVector centerPoint)
	{
		for (final CoordinateVector[] vertices : cell.getSides())
		{
			final boolean isOnBoundary = Arrays.stream(vertices)
			                                   .map(centerPoint::sub)
			                                   .mapToDouble(CoordinateVector::euclidianNorm)
			                                   .allMatch(norm -> DoubleCompare.almostEqual(norm,
			                                                                               radius));
			final DistortedFace face = new DistortedFace(vertices, isOnBoundary);
			for (final DistortedCell c : getCellsBelongingToFace(face, neighbouringCells))
			{
				c.faces.add(face);
				face.addCell(c);
			}
			faces.add(face);
		}
	}
	
	@NotNull
	private static List<DistortedCell> getRefinedNeighbouringCells(final Multimap<DistortedCell, DistortedCell> refinedCells,
	                                                               final DistortedCell cell)
	{
		final List<DistortedCell> neighbouringCells = new ArrayList<>();
		for (final DistortedFace f : cell.faces)
			neighbouringCells.addAll(f.getCells());
		final List<DistortedCell> refinedNeighbouringCells = new ArrayList<>();
		for (final DistortedCell c : neighbouringCells)
			refinedNeighbouringCells.addAll(refinedCells.get(c));
		return refinedNeighbouringCells;
	}
	
	CoordinateVector mean(final CoordinateVector... points)
	{
		return Arrays.stream(points).reduce(new CoordinateVector(dimension), CoordinateVector::add).mul(
			1. / points.length);
	}
	
	List<DistortedCell> partitionCell(final DistortedCell cell)
	{
		final List<DistortedCell> ret = new ArrayList<>();
		final CoordinateVector[] vertices = cell.vertices;
		if (dimension == 2)
		{
			final DistortedFace bottom = cell.getFaceFromVertexNumbers(0, 1);
			CoordinateVector bottomCenter = bottom.center();
			final DistortedFace right = cell.getFaceFromVertexNumbers(1, 2);
			CoordinateVector rightCenter = right.center();
			final DistortedFace top = cell.getFaceFromVertexNumbers(2, 3);
			CoordinateVector topCenter = top.center();
			final DistortedFace left = cell.getFaceFromVertexNumbers(3, 0);
			CoordinateVector leftCenter = left.center();
			if (bottom.isBoundaryFace())
				bottomCenter = centerPoint.add(bottomCenter.sub(centerPoint).normalize().mul(radius));
			if (right.isBoundaryFace())
				rightCenter = centerPoint.add(rightCenter.sub(centerPoint).normalize().mul(radius));
			if (top.isBoundaryFace())
				topCenter = centerPoint.add(topCenter.sub(centerPoint).normalize().mul(radius));
			if (left.isBoundaryFace())
				leftCenter = centerPoint.add(leftCenter.sub(centerPoint).normalize().mul(radius));
			final CoordinateVector center = mean(bottomCenter, rightCenter, topCenter, leftCenter);
			ret.add(new DistortedCell(vertices[0], bottomCenter, center, leftCenter));
			ret.add(new DistortedCell(bottomCenter, vertices[1], rightCenter, center));
			ret.add(new DistortedCell(leftCenter, center, topCenter, vertices[3]));
			ret.add(new DistortedCell(center, rightCenter, vertices[2], topCenter));
			return ret;
		}
		if (dimension == 3)
		{
			throw new UnsupportedOperationException("Not Yet implemented");
		}
		throw new IllegalStateException("Dimension Wrong");
	}
	
	public Multimap<DistortedCell, DistortedCell> refineCells(final Multimap<DistortedCell, DistortedCell> genCells)
	{
		final Multimap<DistortedCell, DistortedCell> refinedCells = HashMultimap.create();
		for (final DistortedCell cell : genCells.values())
			refinedCells.putAll(cell, partitionCell(cell));
		return refinedCells;
	}
	
	public List<CoordinateVector> generatePlotPoints(final int resolution)
	{
		final int pointsPerCell = resolution / cells.size();
		final List<CoordinateVector> ret = new ArrayList<>();
		final IntCoordinates pointsPerDimension = IntCoordinates.repeat(pointsPerCell, dimension);
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
						                coords.set(1. / (pointsPerCell - 1) *
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
