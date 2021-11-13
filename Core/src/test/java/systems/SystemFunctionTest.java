package systems;

import basic.Function;
import basic.FunctionSignature;
import basic.ScalarFunction;
import basic.VectorFunction;
import linalg.*;
import org.junit.Test;

import static org.junit.Assert.assertEquals;

public class SystemFunctionTest
{
	private static SystemFunction createSystemFunction()
	{
		final ScalarFunction f1 = new ScalarFunction()
		{
			@Override
			public int getDomainDimension()
			{
				return 2;
			}
			
			@Override
			public Double value(final CoordinateVector pos)
			{
				return pos.x();
			}
			
			@Override
			public CoordinateVector gradient(final CoordinateVector pos)
			{
				return CoordinateVector.fromValues(1, 0);
			}
		};
		final VectorFunction f2 = new VectorFunction()
		{
			@Override
			public int getRangeDimension()
			{
				return 3;
			}
			
			@Override
			public int getDomainDimension()
			{
				return 2;
			}
			
			@Override
			public CoordinateVector value(final CoordinateVector pos)
			{
				return CoordinateVector.fromValues(pos.x() * pos.x(), pos.y(), 1);
			}
			
			@Override
			public CoordinateMatrix gradient(final CoordinateVector pos)
			{
				return CoordinateDenseMatrix.fromValues(2, 3, 2 * pos.x(), 0, 0, 0, 1, 0);
			}
		};
		final ScalarFunction f3 = new ScalarFunction()
		{
			@Override
			public int getDomainDimension()
			{
				return 2;
			}
			
			@Override
			public Double value(final CoordinateVector pos)
			{
				return pos.y();
			}
			
			@Override
			public CoordinateVector gradient(final CoordinateVector pos)
			{
				return CoordinateVector.fromValues(0, 1);
			}
		};
		final VectorFunction f4 = new VectorFunction()
		{
			@Override
			public int getRangeDimension()
			{
				return 2;
			}
			
			@Override
			public int getDomainDimension()
			{
				return 2;
			}
			
			@Override
			public CoordinateVector value(final CoordinateVector pos)
			{
				return pos;
			}
			
			@Override
			public CoordinateMatrix gradient(final CoordinateVector pos)
			{
				return CoordinateDenseMatrix.fromValues(2, 2, 1, 0, 0, 1);
			}
		};
		return new SystemFunction(new Function[]{f1, f2, f3, f4});
	}
	
	@Test
	public void testValue()
	{
		SystemParameters.createInstance(new int[]{1, 4, 5, 7},
		                                new FunctionSignature[]{new FunctionSignature(Double.class,
		                                                                              CoordinateVector.class,
		                                                                              CoordinateMatrix.class),
		                                                        new FunctionSignature(CoordinateVector.class,
		                                                                              CoordinateMatrix.class,
		                                                                              CoordinateTensor.class),
		                                                        new FunctionSignature(Double.class,
		                                                                              CoordinateVector.class,
		                                                                              CoordinateMatrix.class),
		                                                        new FunctionSignature(CoordinateVector.class,
		                                                                              CoordinateMatrix.class,
		                                                                              CoordinateTensor.class)});
		final SystemFunction f = createSystemFunction();
		assertEquals(f.value(new CoordinateVector(2)), DenseVector.vectorFromValues(0, 0, 0, 1, 0, 0, 0));
		assertEquals(f.value(CoordinateVector.fromValues(1.3, 2)),
		             DenseVector.vectorFromValues(1.3, 1.3 * 1.3, 2, 1, 2, 1.3, 2));
		SystemParameters.deleteInstance();
	}
	
	@Test
	public void testValueSingleComponent()
	{
		SystemParameters.createInstance(new int[]{1, 4, 5, 7},
		                                new FunctionSignature[]{new FunctionSignature(Double.class,
		                                                                              CoordinateVector.class,
		                                                                              CoordinateMatrix.class),
		                                                        new FunctionSignature(CoordinateVector.class,
		                                                                              CoordinateMatrix.class,
		                                                                              CoordinateTensor.class),
		                                                        new FunctionSignature(Double.class,
		                                                                              CoordinateVector.class,
		                                                                              CoordinateMatrix.class),
		                                                        new FunctionSignature(CoordinateVector.class,
		                                                                              CoordinateMatrix.class,
		                                                                              CoordinateTensor.class)});
		final SystemFunction f = new SystemFunction(new ScalarFunction()
		{
			@Override
			public int getDomainDimension()
			{
				return 2;
			}
			
			@Override
			public Double value(final CoordinateVector pos)
			{
				return 7.8;
			}
		}, 2);
		assertEquals(f.value(new CoordinateVector(2)), DenseVector.vectorFromValues(0, 0, 0, 0, 7.8, 0, 0));
		SystemParameters.deleteInstance();
	}
	
	@Test
	public void testGradient()
	{
		SystemParameters.createInstance(new int[]{1, 4, 5, 7},
		                                new FunctionSignature[]{new FunctionSignature(Double.class,
		                                                                              CoordinateVector.class,
		                                                                              CoordinateMatrix.class),
		                                                        new FunctionSignature(CoordinateVector.class,
		                                                                              CoordinateMatrix.class,
		                                                                              CoordinateTensor.class),
		                                                        new FunctionSignature(Double.class,
		                                                                              CoordinateVector.class,
		                                                                              CoordinateMatrix.class),
		                                                        new FunctionSignature(CoordinateVector.class,
		                                                                              CoordinateMatrix.class,
		                                                                              CoordinateTensor.class)});
		final SystemFunction f = createSystemFunction();
		final DenseMatrix d = new DenseMatrix(7, 2);
		d.set(1, 0, 0);
		d.set(0, 0, 1);
		d.set(0, 1, 0);
		d.set(0, 1, 1);
		d.set(0, 2, 0);
		d.set(1, 2, 1);
		d.set(0, 3, 0);
		d.set(0, 3, 1);
		d.set(0, 4, 0);
		d.set(1, 4, 1);
		d.set(1, 5, 0);
		d.set(0, 5, 1);
		d.set(0, 6, 0);
		d.set(1, 6, 1);
		assertEquals(f.gradient(new CoordinateVector(2)), d.transpose());
		SystemParameters.deleteInstance();
	}
}
