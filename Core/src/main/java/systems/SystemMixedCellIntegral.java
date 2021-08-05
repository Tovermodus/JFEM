package systems;

import basic.*;
import linalg.CoordinateMatrix;
import linalg.CoordinateVector;
import tensorproduct.QuadratureRule1D;
import tensorproduct.geometry.TPCell;
import tensorproduct.TPCellIntegral;

public abstract class SystemMixedCellIntegral<CT extends Cell<CT,?>>
	extends CellIntegral<CT, SystemShapeFunction<CT,?, ?>>
{
	protected final int component1;
	protected final int component2;
	
	protected SystemMixedCellIntegral(String name, QuadratureRule1D quadratureRule1D, int component1,
	                                  int component2)
	{
		super(name, quadratureRule1D);
		this.component1 = component1;
		this.component2 = component2;
	}
	
	public SystemMixedCellIntegral(String name, Function<?,?,?> weight, QuadratureRule1D quadratureRule1D, int component1, int component2)
	{
		super(weight, name, quadratureRule1D);
		this.component1 = component1;
		this.component2 = component2;
	}
	
	@Override
	public abstract double evaluateCellIntegral(CT cell, SystemShapeFunction<CT, ?, ?> shapeFunction1,
	                                    SystemShapeFunction<CT, ?, ?> shapeFunction2);
	
}
