package com.simonschmidt;

import java.io.IOException;

public class SparseSolver
{
	static {
		try
		{
			NativeUtils.loadLibraryFromJar("/native/libnative.so");
		} catch (IOException e)
		{
			System.loadLibrary("native");
		}
	}
	final int[] Xs;
	final int[] Ys;
	final double[] Vals;
	final int valueCount;
	final int rows;
	final int cols;
	public SparseSolver(int[] Xs, int [] Ys, double[] Vals, int valueCount, int rows, int cols)
	{
		this.Xs = Xs;
		this.Ys = Ys;
		this.Vals = Vals;
		this.valueCount = valueCount;
		this.rows = rows;
		this.cols = cols;
	}
	public double[] solve(double[] vect)
	{
		return iSolve(Xs, Ys, Vals, valueCount, rows, cols, vect);
	}
	public double[][] inverse()
	{
		double[] inverseFlat = iInverse(Xs, Ys, Vals, valueCount, rows, cols);
		double [][] ret = new double[rows][rows];
		for(int i = 0; i < rows; i++)
			if (cols >= 0) System.arraycopy(inverseFlat, i * rows, ret[i], 0, cols);
		return ret;
	}
	private native double[] iSolve(int[] Xs, int [] Ys, double[] Vals, int valueCount, int rows, int cols,
	                               double[] vect);
	private native double[] iInverse(int[] Xs, int [] Ys, double[] Vals, int valueCount, int rows, int cols);
	public native void sayHello();
	public void hi()
	{
		sayHello();
	}
}
