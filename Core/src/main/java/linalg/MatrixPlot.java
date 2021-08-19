package linalg;

import basic.Plot;

import java.awt.*;

public class MatrixPlot extends Plot
{
	private final Matrix m;
	public int pixelWidth;
	public int pixelHeight;
	public double maxValue;
	final double IDENTIFIED_AS_ZERO = 1e-4;
	
	public MatrixPlot(final Matrix m)
	{
		this.m = m;
		this.maxValue = m.absMaxElement();
		if (maxValue > 1e30)
			maxValue = 1e30;
	}
	
	@Override
	public void drawValues(final Graphics g, final int width, final int height, final double slider)
	{
		pixelWidth = (width - 150) / m.getCols() + 1;
		pixelHeight = (height - 150) / m.getRows() + 1;
		for (final IntCoordinates c : m.getShape().range())
		{
			drawSinglePoint(g, width, height, m.at(c), c);
		}
	}
	
	public void drawSinglePoint(final Graphics g, final int width, final int height, double val, final IntCoordinates coords)
	{
		if (Math.abs(val) > 1e30)
			val = 1e30;
		double logval = 0;
		final double minlog = Math.log(IDENTIFIED_AS_ZERO);
		if (Math.abs(val) > IDENTIFIED_AS_ZERO)
			logval = Math.log(Math.abs(val)) - minlog;
		double maxlog = 1;
		if (maxValue > IDENTIFIED_AS_ZERO)
			maxlog = Math.log(maxValue) - minlog;
		//System.out.println(Math.log(Math.abs(val))+" "+logval+ " " + minlog + " " + maxlog + " " + Math.log
		// (maxValue));
		final int red = (int) (logval / maxlog * 255);
		final int green = 0;
		final int blue = 255 - red;
		final Color c = new Color(red, green, blue);
		final int x = (int) ((1.0 * coords.get(1) / m.getCols()) * (width - 150) + 75 - pixelWidth / 2);
		final int y = (int) ((1.0 * coords.get(0) / m.getRows()) * (height - 150) - pixelWidth / 2) + 75;
		g.setColor(c);
		g.fillRect(x, y, pixelWidth, pixelHeight);
	}
	
	@Override
	public String title()
	{
		return "Matrix. Scale logarithmic: blue = " + IDENTIFIED_AS_ZERO + ", red = " + maxValue;
	}
}
