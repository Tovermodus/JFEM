package systems;

import basic.FunctionSignature;
import linalg.*;
import org.junit.jupiter.api.Test;
import tensorproduct.*;
import tensorproduct.geometry.CartesianGrid;
import tensorproduct.geometry.Cell1D;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

import java.util.List;

import static org.junit.jupiter.api.Assertions.assertEquals;

public class SystemCellIntegralTest
{
	@Test
	public void testMixed() {
		SystemParameters.createInstance(new int[]{1,4,5,7}, new FunctionSignature[]{
			new FunctionSignature(Double.class, CoordinateVector.class, CoordinateMatrix.class),
			new FunctionSignature(CoordinateVector.class, CoordinateMatrix.class, CoordinateTensor.class),
			new FunctionSignature(Double.class, CoordinateVector.class, CoordinateMatrix.class),
			new FunctionSignature(CoordinateVector.class, CoordinateMatrix.class, CoordinateTensor.class)
		});
		CartesianGrid grid = new CartesianGrid(CoordinateVector.fromValues(0,-1),
			CoordinateVector.fromValues(1,2), new IntCoordinates(1,1));
		TPCell c = grid.cells.get(0);
		SystemShapeFunction<TPCell, TPFace, TPShapeFunction> f1
			= new SystemShapeFunction<>(new TPShapeFunction(c,1,0),0);
		SystemShapeFunction<TPCell, TPFace, TPShapeFunction> f2
			= new SystemShapeFunction<>(new TPShapeFunction(c,1,0),2);
		SystemMixedCellIntegral<TPCell> integral =
			new SystemMixedTPCellIntegral(SystemMixedTPCellIntegral.VALUE_VALUE,0,2);
		System.out.println(integral.evaluateCellIntegral(c,f1,f2));
		SystemParameters.deleteInstance();
	}
}
