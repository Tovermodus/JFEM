package basic;

import linalg.MutableMatrix;
import linalg.MutableVector;

import java.util.Set;
import java.util.function.BiConsumer;
import java.util.function.BiPredicate;
import java.util.function.Consumer;
import java.util.function.Predicate;

public interface AcceptsBoundaryValues<CT extends Cell<CT,FT>,
					FT extends Face<CT,FT>,
					ST extends ShapeFunction<CT, FT, valueT, gradientT, hessianT>,
					valueT, gradientT, hessianT>
	extends FESpace<CT,FT,ST>
{
	Set<Integer> getFixedNodeIndices();
	void setFixedNodeValue(int index, double value);
	default void setBoundaryValues(Function<valueT, gradientT, hessianT> boundaryValues, Predicate<FT> faces,
	                               BiPredicate<FT, ST> functions)
	{
		forEachBoundaryFace(F -> {
			if (faces.test(F))
			{
				forEachFunctionOnFace(F,(v) ->
				{
					if (functions.test(F, v) && v.getNodeFunctional().usesFace(F))
					{
						setFixedNodeValue(v.getGlobalIndex(),
							v.getNodeFunctional().evaluate(boundaryValues));
					}
				});
			}
		});
	}
	default void setBoundaryValues(Function<valueT, gradientT, hessianT> boundaryValues, Predicate<FT> faces)
	{
		setBoundaryValues(boundaryValues, faces, (f,st) -> true);
	}
	default void setBoundaryValues(Function<valueT, gradientT, hessianT> boundaryValues)
	{
		setBoundaryValues(boundaryValues, ft -> true, (f,st) -> true);
	}
	@Override
	default void forEachFunctionOnCell(CT K, Consumer<ST> action)
	{
		for (ST function : getShapeFunctionsWithSupportOnCell(K))
		{
			if (getFixedNodeIndices().contains(function.getGlobalIndex()))
				continue;
			action.accept(function);
		}
	}
	
	@Override
	default void forEachFunctionOnFace(FT F, Consumer<ST> action)
	{
		for (ST function : getShapeFunctionsWithSupportOnFace(F))
		{
			if (getFixedNodeIndices().contains(function.getGlobalIndex()))
				continue;
			action.accept(function);
		}
	}
	
	@Override
	default void forEachFunctionCombinationOnCell(CT K, BiConsumer<ST, ST> action)
	{
		for (ST function1 : getShapeFunctionsWithSupportOnCell(K))
		{
			//if (getFixedNodeIndices().contains(function1.getGlobalIndex()))
			//	continue;
			for (ST function2 : getShapeFunctionsWithSupportOnCell(K))
			{
				if (getFixedNodeIndices().contains(function2.getGlobalIndex()))
					continue;
				action.accept(function1, function2);
			}
		}
	}
	
	@Override
	default void forEachFunctionCombinationOnFace(FT F, BiConsumer<ST, ST> action)
	{
		for (ST function1 : getShapeFunctionsWithSupportOnFace(F))
		{
			//if (getFixedNodeIndices().contains(function1.getGlobalIndex()))
			//	continue;
			for (ST function2 : getShapeFunctionsWithSupportOnFace(F))
			{
				if (getFixedNodeIndices().contains(function2.getGlobalIndex()))
					continue;
				action.accept(function1, function2);
			}
		}
	}
	
}
