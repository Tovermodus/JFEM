package basic;

import io.vavr.Function3;
import linalg.MutableMatrix;
import linalg.MutableVector;

import java.util.Set;
import java.util.function.BiConsumer;
import java.util.function.BiFunction;
import java.util.function.Consumer;

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
		getSystemMatrix().deleteRow(index);
		//getSystemMatrix().deleteColumn(index);
		getSystemMatrix().set(1, index, index);
		getRhs().set(value, index);
	}
	
}
