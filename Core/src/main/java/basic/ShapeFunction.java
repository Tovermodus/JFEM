package basic;

import linalg.CoordinateVector;

import java.util.Map;
import java.util.Set;

public interface ShapeFunction<CT extends Cell<CT,FT,ST>, FT extends Face<CT,FT,ST>, ST extends ShapeFunction<CT,FT,ST,
	valueT,gradientT,hessianT>, valueT,	gradientT,
	hessianT> extends Function<valueT,gradientT,hessianT>, Comparable<ST>
{
	@Override
	default int getDomainDimension()
	{
		return getCells().iterator().next().getDimension();
	}
	
	Set<CT> getCells();
	Set<FT> getFaces();
	
	NodeFunctional getNodeFunctional();
	
	void setGlobalIndex(int index);
	int getGlobalIndex();
	
	void addFace(FT face);
	void addCell(CT cell);
	
	valueT valueInCell(CoordinateVector pos, CT cell);
	default gradientT gradientInCell(CoordinateVector pos, CT cell)
	{
		throw new UnsupportedOperationException();
	}
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
	
	Map<Integer, Double> prolongate(Set<ST> refinedFunctions);
	
}
