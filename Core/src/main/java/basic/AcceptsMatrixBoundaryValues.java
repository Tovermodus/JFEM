package basic;

import linalg.MutableMatrix;
import linalg.MutableVector;

import java.util.function.BiPredicate;
import java.util.function.Predicate;

public interface AcceptsMatrixBoundaryValues<CT extends Cell<CT, FT>,
	FT extends Face<CT, FT>,
	ST extends ShapeFunction<CT, FT, valueT, gradientT, hessianT>,
	valueT, gradientT, hessianT>
	extends AcceptsBoundaryValues<CT, FT, ST, valueT, gradientT, hessianT>, MatrixFESpace<CT, FT, ST>
{
	
	@Override
	default void setFixedNodeValue(final int index, final double value)
	{
		getFixedNodeIndices().add(index);
		overWriteValue(index, value, getSystemMatrix(), getRhs());
	}
	
	default void overWriteValue(final int index, final double value, final MutableMatrix mat, final MutableVector rhs)
	{
		mat.deleteRow(index);
		mat.set(1, index, index);
		rhs.set(value, index);
	}
	
	default void writeBoundaryValuesTo(final Function<valueT, gradientT, hessianT> boundaryValues,
	                                   final Predicate<FT> faces,
	                                   final BiPredicate<FT, ST> functions, final MutableMatrix A, final MutableVector rhs)
	{
		forEachBoundaryFace(F ->
		                    {
			                    if (faces.test(F))
			                    {
				                    forEachFunctionOnFace(F, (v) ->
				                    {
					                    if (functions.test(F, v) && v.getNodeFunctional().usesFace(
						                    F))
					                    {
						                    overWriteValue(v.getGlobalIndex(),
						                                   v
							                                   .getNodeFunctional()
							                                   .evaluate(boundaryValues), A,
						                                   rhs);
					                    }
				                    });
			                    }
		                    });
	}
	
	default void writeBoundaryValuesTo(final Function<valueT, gradientT, hessianT> boundaryValues, final Predicate<FT> faces,
	                                   final MutableMatrix A, final MutableVector rhs)
	{
		writeBoundaryValuesTo(boundaryValues, faces, (f, st) -> true, A, rhs);
	}
	
	default void writeBoundaryValuesTo(final Function<valueT, gradientT, hessianT> boundaryValues, final MutableMatrix A,
	                                   final MutableVector rhs)
	{
		writeBoundaryValuesTo(boundaryValues, f -> true, (f, st) -> true, A, rhs);
	}
}
