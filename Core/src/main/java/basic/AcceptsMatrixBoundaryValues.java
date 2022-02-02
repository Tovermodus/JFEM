package basic;

import linalg.DenseVector;
import linalg.MutableMatrix;
import linalg.MutableVector;
import linalg.Vector;

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
	
	default void overWriteValue(final int index,
	                            final double value,
	                            final MutableMatrix mat,
	                            final MutableVector rhs)
	{
		synchronized (this)
		{
			final DenseVector column = mat.getColumn(index);
			for (int i = 0; i < column.getLength(); i++)
			{
				if (getFixedNodeIndices().contains(i))
					continue;
				rhs.add(-column.at(i) * value, i);
			}
			mat.deleteColumn(index);//in rechte seite schreiben
			mat.deleteRow(index);
			mat.set(1, index, index);
			rhs.set(value, index);
		}
	}
	
	default int[] getBoundaryNodes(final Predicate<FT> faces,
	                               final BiPredicate<FT, ST> functions)
	{
		return getFaces().stream()
		                 .filter(Face::isBoundaryFace)
		                 .filter(faces)
		                 .flatMapToInt(F -> getShapeFunctionsWithSupportOnFace(F).stream()
		                                                                         .filter(v -> functions.test(F,
		                                                                                                     v))
		                                                                         .filter(v -> v.getNodeFunctional()
		                                                                                       .usesFace(F))
		                                                                         .mapToInt(ShapeFunction::getGlobalIndex))
		                 .distinct()
		                 .toArray();
	}
	
	default void writeBoundaryValuesTo(final Function<valueT, gradientT, hessianT> boundaryValues,
	                                   final Predicate<FT> faces,
	                                   final BiPredicate<FT, ST> functions,
	                                   final MutableMatrix A,
	                                   final MutableVector rhs)
	{
		final int[] boundaryNodes = getBoundaryNodes(faces, functions);
		for (final int node : boundaryNodes)
			overWriteValue(node,
			               getShapeFunctionMap().get(node)
			                                    .getNodeFunctional()
			                                    .evaluate(boundaryValues), A,
			               rhs);
	}
	
	default void writeBoundaryValuesTo(final Function<valueT, gradientT, hessianT> boundaryValues,
	                                   final Predicate<FT> faces,
	                                   final MutableMatrix A,
	                                   final MutableVector rhs)
	{
		writeBoundaryValuesTo(boundaryValues, faces, (f, st) -> true, A, rhs);
	}
	
	default void writeBoundaryValuesTo(final Function<valueT, gradientT, hessianT> boundaryValues,
	                                   final MutableMatrix A,
	                                   final MutableVector rhs)
	{
		writeBoundaryValuesTo(boundaryValues, f -> true, (f, st) -> true, A, rhs);
	}
	
	default void projectOntoBoundaryValues(final Function<valueT, gradientT, hessianT> boundaryValues,
	                                       final Predicate<FT> faces,
	                                       final BiPredicate<FT, ST> functions,
	                                       final MutableVector vector)
	{
		var shapeFunctionMap = getShapeFunctionMap();
		final int[] boundaryNodes = getBoundaryNodes(faces, functions);
		for (final int node : boundaryNodes)
			vector.set(shapeFunctionMap.get(node)
				                     .getNodeFunctional()
				                     .evaluate(boundaryValues),
				node);
	}
	
	default void copyBoundaryValues(final Vector source, final MutableVector destination,
	                                final Predicate<FT> faces,
	                                final BiPredicate<FT, ST> functions)
	{
		
		final int[] boundaryNodes = getBoundaryNodes(faces, functions);
		for (final int node : boundaryNodes)
			destination.set(source.at(node), node);
	}
	
	default void copyBoundaryValues(final Vector source, final MutableVector destination,
	                                final Predicate<FT> faces)
	{
		copyBoundaryValues(source, destination, faces, (f, st) -> true);
	}
	
	default void copyBoundaryValues(final Vector source, final MutableVector destination)
	{
		copyBoundaryValues(source, destination, f -> true, (f, st) -> true);
	}
	
	default void projectOntoBoundaryValues(final Function<valueT, gradientT, hessianT> boundaryValues,
	                                       final Predicate<FT> faces,
	                                       final MutableVector vector)
	{
		projectOntoBoundaryValues(boundaryValues, faces, (f, st) -> true, vector);
	}
	
	default void projectOntoBoundaryValues(final Function<valueT, gradientT, hessianT> boundaryValues,
	                                       final MutableVector vector)
	{
		projectOntoBoundaryValues(boundaryValues, f -> true, (f, st) -> true, vector);
	}
}
