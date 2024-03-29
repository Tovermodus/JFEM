package linalg;

import basic.Plot;

import java.awt.*;
import java.util.TreeMap;

public class MatrixPlot
	extends Plot
{
	private final Matrix m;
	public int pixelWidth;
	public int pixelHeight;
	public double maxValue;
	final double IDENTIFIED_AS_ZERO = 1e-8;
	private final String title;
	
	public MatrixPlot(final Matrix m)
	{
		this.m = m;
		this.maxValue = m.absMaxElement();
		if (maxValue > 1e30) maxValue = 1e30;
		this.title = "Matrix";
	}
	
	public MatrixPlot(final Matrix m, final String title)
	{
		this.m = m;
		this.maxValue = m.absMaxElement();
		if (maxValue > 1e30) maxValue = 1e30;
		this.title = title;
	}
	
	@Override
	public void drawValues(final Graphics g, final int width, final int height, final double slider)
	{
		if (m.isSparse())
		{
			pixelWidth = (width - 150) / m.getCols() + 2;
			pixelHeight = (height - 150) / m.getRows() + 2;
			g.setColor(Color.BLUE);
			g.fillRect(75 - pixelWidth / 2,
			           75 - pixelHeight / 2,
			           width - 150,
			           height - 150);
			final TreeMap<IntCoordinates, Double> sortedEntries = new TreeMap<>(m.getCoordinateEntryList());
			for (final var v : sortedEntries
				.entrySet())
			{
				drawSinglePoint(g, width, height, v.getValue(), v.getKey());
			}
		} else
		{
			pixelWidth = (width - 150) / Math.min(m.getCols(), 300) + 2;
			pixelHeight = (height - 150) / Math.min(m.getRows(), 300) + 2;
			final int skips = Math.max(1,
			                           m.getShape()
			                            .size() / 90000);
			int i = 0;
			for (final IntCoordinates c : m.getShape()
			                               .range())
			{
				
				if (i++ % skips == 0) drawSinglePoint(g, width, height, m.at(c), c);
			}
		}
	}
	
	public void drawSinglePoint(final Graphics g,
	                            final int width,
	                            final int height,
	                            double val,
	                            final IntCoordinates coords)
	{
		if (Math.abs(val) > 1e30) val = 1e30;
		double logval = 0;
		final double minlog = Math.log(IDENTIFIED_AS_ZERO);
		if (Math.abs(val) > IDENTIFIED_AS_ZERO) logval = Math.log(Math.abs(val)) - minlog;
		double maxlog = 1;
		if (maxValue > IDENTIFIED_AS_ZERO) maxlog = Math.log(maxValue) - minlog;
		//System.out.println(Math.log(Math.abs(val))+" "+logval+ " " + minlog + " " + maxlog + " " + Math.log
		// (maxValue));
		
		final int red;
		final int green;
		if (val >= 0)
		{
			red = (int) (logval / maxlog * 255);
			green = 0;
		} else
		{
			green = (int) (logval / maxlog * 255);
			red = 0;
		}
		final int blue = 255 - red - green;
		final Color c = new Color(red, green, blue);
		final int x = (int) ((1.0 * coords.get(1) / m.getCols()) * (width - 150) + 75 - pixelWidth / 2);
		final int y = (int) ((1.0 * coords.get(0) / m.getRows()) * (height - 150) - pixelHeight / 2) + 75;
		g.setColor(c);
		g.fillRect(x, y, pixelWidth, pixelHeight);
	}
	
	@Override
	public String title()
	{
		return title + ". Scale logarithmic: blue = " + IDENTIFIED_AS_ZERO + ", red = " + maxValue + ", green " +
			"= " + -maxValue;
	}
}
