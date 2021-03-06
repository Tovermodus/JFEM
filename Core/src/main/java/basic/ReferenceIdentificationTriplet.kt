package basic

data class ReferenceCellIdentificationTriplet<CT: CellWithReferenceCell<CT, FT, ET>,
        FT : FaceWithReferenceFace<CT, FT, ET>?,
        ET : EdgeWithReferenceEdge<CT, FT, ET>?,
        ST: ShapeFunctionWithReferenceShapeFunction<CT,FT,ET, *,*,*>>(
    val firstFunction: ST,
    val secondFunction: ST,
    val cell: CT)

data class ReferenceFaceIdentificationTriplet<CT: CellWithReferenceCell<CT, FT, ET>,
        FT : FaceWithReferenceFace<CT, FT, ET>?,
        ET : EdgeWithReferenceEdge<CT, FT, ET>?,
        ST: ShapeFunctionWithReferenceShapeFunction<CT,FT,ET,*, *, *>>(
    val firstFunction: ST,
    val secondFunction: ST,
    val face: FT)

data class ReferenceEdgeIdentificationTriplet<CT: CellWithReferenceCell<CT, FT, ET>,
        FT : FaceWithReferenceFace<CT, FT, ET>?,
        ET : EdgeWithReferenceEdge<CT, FT, ET>?,
        ST: ShapeFunctionWithReferenceShapeFunction<CT,FT,ET, *, *, *>>(
    val firstFunction: ST,
    val secondFunction: ST,
    val edge: ET)
