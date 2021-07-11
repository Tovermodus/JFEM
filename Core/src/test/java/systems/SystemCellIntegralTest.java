package systems;

import basic.Function;
import basic.FunctionSignature;
import basic.ScalarFunction;
import basic.VectorFunction;
import linalg.*;
import org.junit.jupiter.api.Test;
import tensorproduct.*;

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
		Cell1D c1 = new Cell1D(0,1);
		Cell1D c2 = new Cell1D(0,3);
		TPCell c = new TPCell(List.of(c1,c2));
		SystemShapeFunction<TPCell, TPFace, TPEdge, TPShapeFunction> f1
			= new SystemShapeFunction<>(new TPShapeFunction(c,1,0),0);
		SystemShapeFunction<TPCell, TPFace, TPEdge, RTShapeFunction> f2
			= new SystemShapeFunction<>(new RTShapeFunction(c,2,0),3);
		SystemMixedCellIntegral<TPCell> integral = new SystemMixedCellIntegral<>();
		System.out.println(integral.evaluateCellIntegral(c,f1,f2));
		SystemParameters.deleteInstance();
	}
}
