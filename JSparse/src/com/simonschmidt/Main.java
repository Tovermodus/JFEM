package com.simonschmidt;

import java.util.Arrays;

public class Main
{
	public static void main(String[] args) {
		int n = 100;
		int sparseEntries = (int) (n*0.75);
		int[] Xs = new int[n];
		int[] Ys = new int[n];
		double[] Vals = new double[n];
		double[] vect = new double[n];
		int rows = n/10;
		int cols = n/10;
		for(int i = 0; i < sparseEntries; i++)
		{
			if(i < rows)
			{
				Xs[i] = i;
				Ys[i] = i;
			}
			else
			{
				Xs[i] = (int) (Math.random() * cols);
				Ys[i] = (int) (Math.random() * rows);
			}
			Vals[i] = 1+0.1*Math.random();
			vect[i] = 0.1*Math.random();
		}
		SparseSolver s = new SparseSolver(Xs, Ys, Vals, sparseEntries, rows, cols);
		System.out.println(s.solve(vect));
		System.out.println(Arrays.toString(s.solve(vect)));
		System.out.println(Arrays.deepToString(s.inverse()));
	}
}
