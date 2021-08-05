package basic;

import com.google.common.collect.Lists;
import com.google.common.collect.Multimap;
import io.vavr.Function3;
import linalg.CoordinateVector;
import linalg.MutableVector;

import java.util.*;
import java.util.function.BiConsumer;
import java.util.function.BiFunction;
import java.util.function.Consumer;
import java.util.stream.Stream;

public interface FESpace<CT extends Cell<CT,FT>, FT extends  Face<CT,FT>,
	ST extends ShapeFunction<CT,FT,?,?,?>>
{
	int getDimension();
	List<CT> getCells();
	Map<Integer, ST> getShapeFunctions();
	List<FT> getFaces();
	default List<FT> getBoundaryFaces()
	{
		List<FT> boundaryFaces = new ArrayList<>();
		for (FT F : getFaces())
			if (F.isBoundaryFace())
				boundaryFaces.add(F);
		return boundaryFaces;
	}
	Collection<ST> getShapeFunctionsWithSupportOnCell(CT cell);
	Collection<ST> getShapeFunctionsWithSupportOnFace(FT face);
	
	List<CoordinateVector> generatePlotPoints(int resolution);
	default FESpace<CT,FT, ST> refine(Multimap<CT, CT> cellRefinedCellMapping,
	           Multimap<FT, FT> faceRefinedFaceMapping)
	{
		throw new UnsupportedOperationException();
	}
	default void forEachFace(Consumer<FT> action){
		List<List<FT>> smallerList = Lists.partition(getFaces(),getFaces().size()/PerformanceArguments.getInstance().threadNumber+1);
		Stream<List<FT>> stream = smallerList.stream();
		if(PerformanceArguments.getInstance().parallelizeThreads)
			stream = stream.parallel();
		stream.forEach(smallList->
		{
			for (FT F : smallList)
			{
				action.accept(F);
			}
		});
	}
	default void forEachBoundaryFace(Consumer<FT> action){
		forEachFace(F->{
			if(F.isBoundaryFace())
				action.accept(F);
		});
	}
	default void forEachCell(Consumer<CT> action){
		List<List<CT>> smallerList = Lists.partition(getCells(),
			getFaces().size()/PerformanceArguments.getInstance().threadNumber+1);
		Stream<List<CT>> stream = smallerList.stream();
		if(PerformanceArguments.getInstance().parallelizeThreads)
			stream = stream.parallel();
		stream.forEach(smallList->
		{
			for (CT K: smallList)
			{
				action.accept(K);
			}
		});
	}
	default void forEachFunctionOnCell(CT K, Consumer<ST> action){
		for(ST function: getShapeFunctionsWithSupportOnCell(K))
			action.accept(function);
	}
	default void forEachFunctionOnFace(FT F, Consumer<ST> action){
		for(ST function: getShapeFunctionsWithSupportOnFace(F))
			action.accept(function);
	}
	default void forEachFunctionCombinationOnCell(CT K, BiConsumer<ST,ST> action){
		for(ST function1: getShapeFunctionsWithSupportOnCell(K))
			for(ST function2: getShapeFunctionsWithSupportOnCell(K))
				action.accept(function1, function2);
	}
	default void forEachFunctionCombinationOnFace(FT F, BiConsumer<ST,ST> action){
		for(ST function1: getShapeFunctionsWithSupportOnFace(F))
			for(ST function2: getShapeFunctionsWithSupportOnFace(F))
				action.accept(function1, function2);
	}
	
}
