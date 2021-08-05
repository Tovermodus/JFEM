package basic;

public interface FaceWithReferenceFace<CT extends CellWithReferenceCell<CT, FT>, FT extends FaceWithReferenceFace<CT, FT>
	> extends Face<CT,FT>
{
	FT getReferenceFace();
}
