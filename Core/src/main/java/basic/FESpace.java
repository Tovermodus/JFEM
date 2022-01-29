package basic;

import com.google.common.collect.Multimap;
import linalg.CoordinateVector;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.function.BiConsumer;
import java.util.function.Consumer;
import java.util.stream.Stream;

public interface FESpace<CT extends Cell<CT, FT>, FT extends Face<CT, FT>,
	ST extends ShapeFunction<CT, FT, ?, ?, ?>>
{
	int getDimension();
	
	List<CT> getCells();
	
	Map<Integer, ST> getShapeFunctionMap();
	
	Collection<ST> getShapeFunctions();
	
	List<FT> getFaces();
	
	default List<FT> getBoundaryFaces()
	{
		final List<FT> boundaryFaces = new ArrayList<>();
		for (final FT F : getFaces())
			if (F.isBoundaryFace())
				boundaryFaces.add(F);
		return boundaryFaces;
	}
	
	Collection<ST> getShapeFunctionsWithSupportOnCell(CT cell);
	
	Collection<ST> getShapeFunctionsWithSupportOnFace(FT face);
	
	List<CoordinateVector> generatePlotPoints(int resolution);
	
	default FESpace<CT, FT, ST> refine(final Multimap<CT, CT> cellRefinedCellMapping,
	                                   final Multimap<FT, FT> faceRefinedFaceMapping)
	{
		throw new UnsupportedOperationException();
	}
	
	default void forEachFace(final Consumer<FT> action)
	{
		Stream<FT> stream = getFaces().stream();
		if (PerformanceArguments.getInstance().parallelizeThreads)
			stream = stream.parallel();
		stream.forEach(action);
	}
	
	default void forEachBoundaryFace(final Consumer<FT> action)
	{
		forEachFace(F ->
		            {
			            if (F.isBoundaryFace())
				            action.accept(F);
		            });
	}
	
	default void forEachCell(final Consumer<CT> action)
	{
		Stream<CT> stream = getCells().stream();
		if (PerformanceArguments.getInstance().parallelizeThreads)
			stream = stream.parallel();
		stream.forEach(action);
	}
	
	default void forEachFunctionOnCell(final CT K, final Consumer<ST> action)
	{
		for (final ST function : getShapeFunctionsWithSupportOnCell(K))
			action.accept(function);
	}
	
	default void forEachFunctionOnFace(final FT F, final Consumer<ST> action)
	{
		for (final ST function : getShapeFunctionsWithSupportOnFace(F))
			action.accept(function);
	}
	
	default void forEachFunctionCombinationOnCell(final CT K, final BiConsumer<ST, ST> action)
	{
		for (final ST function1 : getShapeFunctionsWithSupportOnCell(K))
			for (final ST function2 : getShapeFunctionsWithSupportOnCell(K))
			{
				action.accept(function1, function2);
			}
	}
	
	default void forEachFunctionCombinationOnFace(final FT F, final BiConsumer<ST, ST> action)
	{
		for (final ST function1 : getShapeFunctionsWithSupportOnFace(F))
			for (final ST function2 : getShapeFunctionsWithSupportOnFace(F))
				action.accept(function1, function2);
	}
	
	Multimap<CT, ST> getCellSupportMapping();
	
	Multimap<FT, ST> getFaceSupportMapping();
}
