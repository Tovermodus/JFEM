package basic

data class ReferenceCellIdentificationTriplet<CT: CellWithReferenceCell<CT, FT>,
        FT : FaceWithReferenceFace<CT, FT>?,
        ST: ShapeFunctionWithReferenceShapeFunction<CT,FT, *,*,*>>(
    val firstFunction: ST,
    val secondFunction: ST,
    val cell: CT)

data class ReferenceFaceIdentificationTriplet<CT: CellWithReferenceCell<CT, FT>,
        FT : FaceWithReferenceFace<CT, FT>?,
        ST: ShapeFunctionWithReferenceShapeFunction<CT,FT,*, *, *>>(
    val firstFunction: ST,
    val secondFunction: ST,
    val face: FT)

data class ReferenceEdgeIdentificationTriplet<CT: CellWithReferenceCell<CT, FT>,
        FT : FaceWithReferenceFace<CT, FT>?,
        ST: ShapeFunctionWithReferenceShapeFunction<CT,FT, *, *, *>>(
    val firstFunction: ST,
    val secondFunction: ST)
