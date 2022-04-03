package tensorproduct.geometry;

import basic.DoubleCompare;
import basic.FaceWithReferenceFace;
import basic.PerformanceArguments;
import basic.VectorFunction;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableSet;
import com.google.common.collect.Multimap;
import linalg.CoordinateComparator;
import linalg.CoordinateVector;
import org.jetbrains.annotations.NotNull;

import java.util.*;

public class TPFace
	implements FaceWithReferenceFace<TPCell, TPFace>, Comparable<TPFace>
{
	public final double otherCoordinate;
	public final int flatDimension;
	final ImmutableList<Cell1D> cell1Ds;
	private final boolean isBoundaryFace;
	private final transient VectorFunction normal;
	final Set<TPCell> cells;
	
	TPFace(final List<Cell1D> cell1Ds,
	       final int flatDimension,
	       final double otherCoordinate,
	       final boolean isBoundaryFace)
	{
		this.otherCoordinate = otherCoordinate;
		this.flatDimension = flatDimension;
		this.isBoundaryFace = isBoundaryFace;
		this.cell1Ds = ImmutableList.copyOf(cell1Ds);
		this.cells = new HashSet<>();
		this.normal = new VectorFunction()
		{
			@Override
			public int getRangeDimension()
			{
				return getDimension();
			}
			
			@Override
			public int getDomainDimension()
			{
				return getDimension();
			}
			
			@Override
			public CoordinateVector value(final CoordinateVector pos)
			{
				final CoordinateVector ret = new CoordinateVector(pos.getLength());
				ret.set(1, flatDimension);
				return ret;
			}
		};
	}
	
	public static TPFace fromVertices(final CoordinateVector[] vertices, final boolean isBoundaryFace)
	{
		final int dimension = vertices[0].getLength();
		int flatDimension = 0;
		double otherCoordinate = 0;
		final List<Cell1D> cell1DS = new ArrayList<>();
		for (int i = 0; i < dimension; i++)
		{
			final int finalI = i;
			final double min = Arrays.stream(vertices)
			                         .mapToDouble(v -> v.at(finalI))
			                         .min()
			                         .getAsDouble();
			final double max = Arrays.stream(vertices)
			                         .mapToDouble(v -> v.at(finalI))
			                         .max()
			                         .getAsDouble();
			if (DoubleCompare.almostEqual(min, max))
			{
				flatDimension = i;
				otherCoordinate = min;
			} else
			{
				cell1DS.add(new Cell1D(min, max));
			}
		}
		return new TPFace(cell1DS, flatDimension, otherCoordinate, isBoundaryFace);
	}
	
	public static TPFace unitHyperCubeFace(final int dimension, final boolean isBoundaryFace)
	{
		if (dimension == 2)
		{
			return new TPFace(List.of(new Cell1D(0, 1)), 1, 0, isBoundaryFace);
		}
		final List<Cell1D> cells = new ArrayList<>(dimension - 1);
		for (int i = 0; i < dimension - 1; i++)
			cells.add(new Cell1D(0, 1));
		return new TPFace(cells, 0, 0, isBoundaryFace);
	}
	
	public ImmutableList<Cell1D> getComponentCells()
	{
		return cell1Ds;
	}
	
	@Override
	public int getDimension()
	{
		return cell1Ds.size() + 1;
	}
	
	@Override
	public ImmutableSet<TPCell> getCells()
	{
		return ImmutableSet.copyOf(cells);
	}
	
	@Override
	public boolean isBoundaryFace()
	{
		return isBoundaryFace;
	}
	
	@Override
	public VectorFunction getNormal()
	{
		return normal;
	}
	
	public boolean isNormalDownstream(final CoordinateVector pos)
	{
		return getNormal().value(center())
		                  .inner(pos.sub(center())) > 0;
	}
	
	public TPCell getUpStreamCell(final CoordinateVector direction)
	{
		if (normal.value(center())
		          .inner(direction) > 0)
			return getNormalUpstreamCell();
		else
			return getNormalDownstreamCell();
	}
	
	public TPCell getDownStreamCell(final CoordinateVector direction)
	{
		if (normal.value(center())
		          .inner(direction) < 0)
			return getNormalUpstreamCell();
		else
			return getNormalDownstreamCell();
	}
	
	@Override
	public TPCell getNormalDownstreamCell()
	{
		for (final TPCell cell : cells)
		{
			if (cell.center()
			        .at(flatDimension) > otherCoordinate)
			{
				return cell;
			}
		}
		return null;
	}
	
	@Override
	public TPCell getNormalUpstreamCell()
	{
		for (final TPCell cell : cells)
		{
			if (cell.center()
			        .at(flatDimension) < otherCoordinate)
			{
				return cell;
			}
		}
		return null;
	}
	
	@Override
	public CoordinateVector center()
	{
		final CoordinateVector ret = new CoordinateVector(getDimension());
		int subd = 0;
		for (int d = 0; d < getDimension(); d++)
		{
			if (d == flatDimension)
				ret.set(otherCoordinate, d);
			else
				ret.set(cell1Ds.get(subd++)
				               .center(), d);
		}
		return ret;
	}
	
	@Override
	public boolean isOnFace(final CoordinateVector pos)
	{
		int subd = 0;
		for (int d = 0; d < getDimension(); d++)
		{
			if (d == flatDimension)
			{
				if (!DoubleCompare.almostEqual(pos.at(flatDimension), otherCoordinate, 1000))
					return false;
			} else
			{
				if (!cell1Ds.get(subd++)
				            .isInCell(pos.at(d)))
					return false;
			}
		}
		return true;
	}
	
	@Override
	public List<TPFace> refine(final Multimap<TPCell, TPCell> cellMap)
	{
		throw new UnsupportedOperationException();
	}
	
	@Override
	public String toString()
	{
		String ret = "";
		int subd = 0;
		for (int d = 0; d < getDimension(); d++)
		{
			if (d == flatDimension)
				ret = ret.concat(otherCoordinate + "");
			else
				ret = ret.concat(
					"[" + cell1Ds.get(subd)
					             .getStart() + ", " + cell1Ds.get(subd++)
					                                         .getEnd() +
						"]");
			if (d < cell1Ds.size())
				ret = ret.concat("x");
		}
		if (isBoundaryFace)
			ret = ret.concat(" On Boundary");
		return ret;
	}
	
	@Override
	public int compareTo(@NotNull final TPFace o)
	{
		if (o.getDimension() < getDimension())
			return -1;
		if (o.getDimension() > getDimension())
			return 1;
		if (o.flatDimension < flatDimension)
			return -1;
		if (o.flatDimension > flatDimension)
			return 1;
		if (o.otherCoordinate < otherCoordinate - PerformanceArguments.getInstance().doubleTolerance)
			return -1;
		if (o.otherCoordinate > otherCoordinate + PerformanceArguments.getInstance().doubleTolerance)
			return 1;
		if (o.isBoundaryFace && !isBoundaryFace)
			return 1;
		if (!o.isBoundaryFace && isBoundaryFace)
			return -1;
		if (o.getNormalDownstreamCell() != null && getNormalDownstreamCell() == null)
			return 1;
		if (o.getNormalDownstreamCell() != null && getNormalDownstreamCell() == null)
			return -1;
		if (o.getNormalUpstreamCell() != null && getNormalUpstreamCell() == null)
			return 1;
		if (o.getNormalUpstreamCell() != null && getNormalUpstreamCell() == null)
			return -1;
		return CoordinateComparator.comp(center().getEntries(),
		                                 o.center()
		                                  .getEntries());
	}
	
	@Override
	public int hashCode()
	{
		int ret = 0;
		for (int i = 0; i < cell1Ds.size(); i++)
		{
			ret += Math.pow(7, 3 * i) * DoubleCompare.doubleHash(cell1Ds.get(i)
			                                                            .center());
			ret += Math.pow(7, 3 * i + 1) * DoubleCompare.doubleHash(cell1Ds.get(i)
			                                                                .getStart());
			ret += Math.pow(7, 3 * i + 2) * DoubleCompare.doubleHash(cell1Ds.get(i)
			                                                                .getEnd());
		}
		ret -= 141 * DoubleCompare.doubleHash(otherCoordinate);
		ret *= (flatDimension + 19);
		ret *= 2;
		if (isBoundaryFace)
			ret += 1;
		if (getNormalDownstreamCell() != null)
			ret += 1789821;
		if (getNormalUpstreamCell() != null)
			ret += 31789821;
		return ret;
	}
	
	@Override
	public boolean equals(final Object obj)
	{
		if (obj instanceof TPFace)
			return 0 == compareTo((TPFace) obj);
		else
			return false;
	}
	
	@Override
	public TPFace getReferenceFace()
	{
		List<Cell1D> cells1 = new ArrayList<>(getDimension() - 1);
		for (int i = 0; i < getDimension() - 1; i++)
		{
			cells1.add(new Cell1D(0, 1));
		}
		final TPFace refFace =
			new TPFace(cells1, flatDimension, 0, isBoundaryFace);
		cells1 = new ArrayList<>(getDimension() - 1);
		final List<Cell1D> cells2 = new ArrayList<>(getDimension() - 1);
		for (int i = 0; i < getDimension(); i++)
		{
			if (i != flatDimension)
			{
				cells1.add(new Cell1D(0, 1));
				cells2.add(new Cell1D(0, 1));
			} else
			{
				cells1.add(new Cell1D(0, 1));
				cells2.add(new Cell1D(-1, 0));
			}
		}
		final TPCell downstreamCell = new TPCell(cells1);
		if (getNormalDownstreamCell() != null)
			refFace.cells.add(downstreamCell);
		final TPCell upstreamCell = new TPCell(cells2);
		if (getNormalUpstreamCell() != null)
			refFace.cells.add(upstreamCell);
		return refFace;
	}
}
