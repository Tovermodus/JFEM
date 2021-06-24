package basic;

public interface FaceWithReferenceFace<CT extends CellWithReferenceCell<CT, FT, ET>, FT extends FaceWithReferenceFace<CT, FT, ET>,
	ET extends EdgeWithReferenceEdge<CT, FT, ET>> extends Face<CT,FT,ET>
{
	public FT getReferenceFace();
}
