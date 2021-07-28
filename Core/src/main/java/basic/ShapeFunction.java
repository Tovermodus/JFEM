package basic;

import com.google.common.collect.Iterables;
import com.google.common.collect.Sets;
import linalg.CoordinateVector;
import org.jetbrains.annotations.NotNull;

import java.util.Map;
import java.util.Set;

public interface ShapeFunction<CT extends Cell<CT,FT, ET>, FT extends Face<CT,FT, ET>,
	ET extends Edge<CT,FT,ET>, valueT,	gradientT,
	hessianT> extends Function<valueT,gradientT,hessianT>
{
	@Override
	default int getDomainDimension()
	{
		return Iterables.getLast(getCells()).getDimension();
	}
	
	Set<CT> getCells();
	
	Set<FT> getFaces();
	
	NodeFunctional<valueT, gradientT, hessianT> getNodeFunctional();
	
	int getGlobalIndex();
	
	valueT valueInCell(CoordinateVector pos, CT cell);
	
	gradientT gradientInCell(CoordinateVector pos, CT cell);
	
	default hessianT hessianInCell(CoordinateVector pos, CT cell)
	{
		throw new UnsupportedOperationException();
	}
	
	valueT jumpInValue(FT face, CoordinateVector pos);
	
	gradientT jumpInDerivative(FT face, CoordinateVector pos);
	
	valueT averageInValue(FT face, CoordinateVector pos);
	
	gradientT averageInDerivative(FT face, CoordinateVector pos);
	
	gradientT normalAverageInValue(FT face, CoordinateVector pos);
	
	valueT normalAverageInDerivative(FT face, CoordinateVector pos);
	
	<ST extends ShapeFunction<CT,FT,ET,valueT,gradientT,hessianT>> Map<Integer, Double> prolongate(Set<ST> refinedFunctions);
	
	int maxPolynomialDegree();
}
