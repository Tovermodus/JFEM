package linalg;

import basic.Plot;
import basic.ScalarFunction;
import com.google.common.collect.Iterables;

import java.awt.*;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.Optional;

public class MatrixPlot implements Plot
{
	private final Matrix m;
	public int pixelWidth;
	public int pixelHeight;
	public double maxValue;
	final double IDENTIFIED_AS_ZERO = 1e-4;
	public MatrixPlot(Matrix m)
	{
		this.m = m;
		this.maxValue = m.absMaxElement();
		if(maxValue > 1e30)
			maxValue = 1e30;
	}
	@Override
	public void drawValues(Graphics g, int width, int height, double slider)
	{
		pixelWidth = (width-150)/m.getCols()+1;
		pixelHeight = (height-150)/m.getRows()+1;
		for(IntCoordinates c: m.getShape().range())
		{
			drawSinglePoint(g, width, height, m.at(c), c);
		}
	}
	public void drawSinglePoint(Graphics g, int width, int height, double val, IntCoordinates coords)
	{
		if(Math.abs(val) > 1e30)
			val = 1e30;
		double logval = 0;
		double minlog = Math.log(IDENTIFIED_AS_ZERO);
		if(Math.abs(val) > IDENTIFIED_AS_ZERO)
			logval = Math.log(Math.abs(val)) - minlog;
		double maxlog = 1;
		if(maxValue > IDENTIFIED_AS_ZERO)
			maxlog = Math.log(maxValue) - minlog;
		//System.out.println(Math.log(Math.abs(val))+" "+logval+ " " + minlog + " " + maxlog + " " + Math.log
		// (maxValue));
		int red =  (int)(logval/maxlog*255);
		int green = 0;
		int blue = 255 - red;
		Color c = new Color(red,green,blue);
		int x = (int)((1.0*coords.get(1)/m.getCols())*(width -150)+75 - pixelWidth/2);
		int y = (int)((1.0*coords.get(0)/m.getRows())*(height -150) - pixelWidth/2) + 75;
		g.setColor(c);
		g.fillRect(x,y,pixelWidth,pixelHeight);
	}
	
	@Override
	public String title()
	{
		return "Matrix. Scale logarithmic: blue = "+IDENTIFIED_AS_ZERO+", red = "+ maxValue;
	}
}
