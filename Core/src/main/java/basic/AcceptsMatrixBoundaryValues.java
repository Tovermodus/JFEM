package basic;

import linalg.MutableMatrix;
import linalg.MutableVector;

import java.util.function.*;

public interface AcceptsMatrixBoundaryValues<CT extends Cell<CT, FT>,
	FT extends Face<CT, FT>,
	ST extends ShapeFunction<CT, FT, valueT, gradientT, hessianT>,
	valueT, gradientT, hessianT>
	extends AcceptsBoundaryValues<CT, FT, ST, valueT, gradientT, hessianT>, MatrixFESpace<CT,FT,ST>
{
	
	@Override
	default void setFixedNodeValue(int index, double value)
	{
		getFixedNodeIndices().add(index);
		overWriteValue(index, value, getSystemMatrix(), getRhs());
	}
	
	default void overWriteValue(int index, double value, MutableMatrix mat, MutableVector rhs)
	{
		mat.deleteRow(index);
		mat.set(1, index, index);
		rhs.set(value, index);
	}
	
	default void writeBoundaryValuesTo(Function<valueT, gradientT, hessianT> boundaryValues,
	                                       Predicate<FT> faces,
	                                 BiPredicate<FT, ST> functions, MutableMatrix A, MutableVector rhs)
	{
		forEachBoundaryFace(F ->
		{
			if (faces.test(F))
			{
				forEachFunctionOnFace(F, (v) ->
				{
					if (functions.test(F, v) && v.getNodeFunctional().usesFace(F))
					{
						overWriteValue(v.getGlobalIndex(),
							v.getNodeFunctional().evaluate(boundaryValues), A, rhs);
					}
				});
			}
		});
	}
	
	default void writeBoundaryValuesTo(Function<valueT, gradientT, hessianT> boundaryValues, Predicate<FT> faces,
	                               MutableMatrix A, MutableVector rhs)
	{
		writeBoundaryValuesTo(boundaryValues, faces, (f, st)->true, A, rhs);
	}
	
	default void writeBoundaryValuesTo(Function<valueT, gradientT, hessianT> boundaryValues, MutableMatrix A,
	                               MutableVector rhs)
	{
		writeBoundaryValuesTo(boundaryValues, f->true, (f, st)->true, A, rhs);
	}
	
}
