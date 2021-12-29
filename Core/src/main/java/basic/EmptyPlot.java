package basic;

import java.awt.*;

public class EmptyPlot
	extends Plot
{
	@Override
	public void drawValues(final Graphics g, final int width, final int height, final double slider)
	{
	
	}
	
	@Override
	public String title()
	{
		return "empty plot";
	}
}
