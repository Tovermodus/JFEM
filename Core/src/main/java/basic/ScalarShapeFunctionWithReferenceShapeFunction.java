package basic;

import linalg.CoordinateMatrix;
import linalg.CoordinateVector;
import org.jetbrains.annotations.NotNull;

import java.util.Set;

public interface ScalarShapeFunctionWithReferenceShapeFunction<CT extends CellWithReferenceCell<CT,FT,ET>,
	FT extends FaceWithReferenceFace<CT,FT,ET>, ET extends EdgeWithReferenceEdge<CT,FT,ET>,
	ST extends ScalarShapeFunctionWithReferenceShapeFunction<CT,FT,ET,ST>> extends ShapeFunctionWithReferenceShapeFunction<CT
	,FT,ET,ST,Double, CoordinateVector, CoordinateMatrix>, ScalarShapeFunction<CT,FT,ET,ST>
{
}
