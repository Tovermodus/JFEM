package distorted.geometry;

import basic.Face;
import basic.FaceWithReferenceFace;
import basic.VectorFunction;
import com.google.common.collect.ImmutableSet;
import com.google.common.collect.Multimap;
import linalg.CoordinateMatrix;
import linalg.CoordinateVector;
import org.jetbrains.annotations.NotNull;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

import java.util.List;
import java.util.Set;

public class DistortedFace implements FaceWithReferenceFace<DistortedCell, DistortedFace>
{
	Set<DistortedCell> cells;
	private final boolean isBoundaryFace;
	private final VectorFunction normal;
	final DistortedFace referenceFace;
	final int dimension;
	
	public DistortedFace()
	{
		isBoundaryFace = false;
		normal = null;
		referenceFace = null;
		dimension = 0;
	}
	
	@Override
	public int getDimension()
	{
		return dimension;
	}
	
	@Override
	public ImmutableSet<DistortedCell> getCells()
	{
		return ImmutableSet.copyOf(cells);
	}
	public CoordinateMatrix getRotationOfDownStreamCell()
	{
		throw new UnsupportedOperationException("not implemented yet");
	
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
	
	@Override
	public DistortedCell getNormalDownstreamCell()
	{
		throw new UnsupportedOperationException("not implemented yet");
	}
	
	@Override
	public DistortedCell getNormalUpstreamCell()
	{
		throw new UnsupportedOperationException("not implemented yet");
	}
	
	@Override
	public CoordinateVector center()
	{
		return getNormalUpstreamCell().transformFromReferenceCell(referenceFace.center());
	}
	
	@Override
	public boolean isOnFace(CoordinateVector pos)
	{
		
		return referenceFace.isOnFace(getNormalUpstreamCell().transformToReferenceCell(pos));
	}
	
	@Override
	public List<DistortedFace> refine(Multimap<DistortedCell, DistortedCell> cellMap)
	{
		throw new UnsupportedOperationException("not implemented yet");
	}
	
	@Override
	public int compareTo(@NotNull DistortedFace o)
	{
		throw new UnsupportedOperationException("not implemented yet");
	}
	
	@Override
	public DistortedFace getReferenceFace()
	{
		return referenceFace;
	}
}
