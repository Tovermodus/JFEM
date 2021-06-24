package basic;

public interface EdgeWithReferenceEdge<CT extends CellWithReferenceCell<CT, FT, ET>, FT extends FaceWithReferenceFace<CT, FT, ET>,
	ET extends EdgeWithReferenceEdge<CT, FT, ET>> extends Edge<CT,FT,ET>
{
	public ET getReferenceEdge();
}
