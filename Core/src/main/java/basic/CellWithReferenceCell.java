package basic;

public interface CellWithReferenceCell<CT extends CellWithReferenceCell<CT, FT, ET>, FT extends FaceWithReferenceFace<CT, FT, ET>,
	ET extends EdgeWithReferenceEdge<CT, FT, ET>> extends Cell<CT,FT,ET>
{
	CT getReferenceCell();
	
	double getGradGradConversionFactor();
}
