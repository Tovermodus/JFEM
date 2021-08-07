package distorted.geometry;

import basic.Face;
import basic.FaceWithReferenceFace;
import basic.VectorFunction;
import com.google.common.collect.ImmutableSet;
import com.google.common.collect.Multimap;
import linalg.CoordinateVector;
import org.jetbrains.annotations.NotNull;

import java.util.List;

public class DistortedFace implements FaceWithReferenceFace<DistortedCell, DistortedFace>
{
	@Override
	public int getDimension()
	{
		throw new UnsupportedOperationException("not implemented yet");
	}
	
	@Override
	public ImmutableSet<DistortedCell> getCells()
	{
		throw new UnsupportedOperationException("not implemented yet");
	}
	
	@Override
	public boolean isBoundaryFace()
	{
		throw new UnsupportedOperationException("not implemented yet");
	}
	
	@Override
	public VectorFunction getNormal()
	{
		throw new UnsupportedOperationException("not implemented yet");
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
		throw new UnsupportedOperationException("not implemented yet");
	}
	
	@Override
	public boolean isOnFace(CoordinateVector pos)
	{
		throw new UnsupportedOperationException("not implemented yet");
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
		throw new UnsupportedOperationException("not implemented yet");
	}
}
