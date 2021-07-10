package systems;

import basic.FunctionSignature;
import linalg.CoordinateMatrix;
import linalg.CoordinateTensor;
import linalg.CoordinateVector;

public class SystemParameters
{
	private static SystemParameters INSTANCE;
	final int[] ends;
	final int[] starts;
	public final FunctionSignature[] signatures;
	
	private SystemParameters(int[] ends, FunctionSignature[] signatures)
	{
		starts = new int[ends.length];
		this.ends = ends.clone();
		this.signatures = signatures;
		for (int i = 0; i < ends.length - 1; i++)
		{
			starts[i + 1] = ends[i];
		}
		for (int i = 0; i < ends.length; i++)
		{
			if(!Double.class.isAssignableFrom(signatures[i].getValueT()) &&  !CoordinateVector.class.isAssignableFrom(signatures[i].getValueT()))
				throw new UnsupportedOperationException("Only Double and CoordinateVector implemented" +
					" as value");
			if(!CoordinateVector.class.isAssignableFrom(signatures[i].getGradientT()) &&  ! CoordinateMatrix.class.isAssignableFrom(signatures[i].getGradientT()))
				throw new UnsupportedOperationException("Only CoordinateVector and CoordinateMatrix " +
					"implemented" +
					" as gradient");
			if(!CoordinateMatrix.class.isAssignableFrom(signatures[i].getHessianT()) &&  ! CoordinateTensor.class.isAssignableFrom(signatures[i].getHessianT()))
				throw new UnsupportedOperationException("Only CoordinateMatrix and CoordinateTensor implemented" +
					" as hessian");
		}
	}
	public static void createInstance(int[] ends, FunctionSignature[] signatures)
	{
		if (INSTANCE == null)
			INSTANCE = new SystemParameters(ends, signatures);
		else
			throw new IllegalStateException("Parameters already set");
	}
	public static void createInstance(int[] ends)
	{
		FunctionSignature[] signats = new FunctionSignature[ends.length];
		for (int i = 0; i < ends.length; i++)
		{
			signats[i] =new FunctionSignature(Double.class,
				CoordinateVector.class,
				CoordinateMatrix.class);
		}
		if (INSTANCE == null)
			INSTANCE = new SystemParameters(ends, signats);
		else
			throw new IllegalStateException("Parameters already set");
	}
	
	public static void deleteInstance()
	{
		INSTANCE = null;
	}
	
	public static SystemParameters getInstance()
	{
		if (INSTANCE == null)
			throw new IllegalStateException("No Parameters set");
		return INSTANCE;
	}
	
	
}
