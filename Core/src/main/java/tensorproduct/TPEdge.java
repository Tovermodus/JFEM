package tensorproduct;

import basic.Edge;
import basic.VectorFunction;
import linalg.CoordinateVector;
import org.jetbrains.annotations.NotNull;

import java.util.Set;

public class TPEdge implements Edge<TPCell, TPFace, TPEdge>
{
	
	@Override
	public int getDimension()
	{
		throw new UnsupportedOperationException("not implemented yet");
	}
	
	@Override
	public Set<TPCell> getCells()
	{
		throw new UnsupportedOperationException("not implemented yet");
	}
	
	@Override
	public VectorFunction getTangent()
	{
		throw new UnsupportedOperationException("not implemented yet");
	}
	
	@Override
	public CoordinateVector center()
	{
		throw new UnsupportedOperationException("not implemented yet");
	}
	
	@Override
	public boolean isOnEdge(CoordinateVector pos)
	{
		throw new UnsupportedOperationException("not implemented yet");
	}
	
	@Override
	public int compareTo(@NotNull TPEdge o)
	{
		throw new UnsupportedOperationException("not implemented yet");
	}
}
