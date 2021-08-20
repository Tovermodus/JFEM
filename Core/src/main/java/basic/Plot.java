package basic;

import java.awt.*;

public abstract class Plot
{
	Overlay o;
	
	void draw(final Graphics g, final int width, final int height, final double slider, final boolean overlay)
	{
		g.setColor(Color.BLACK);
		drawValues(g, width, height - 50, slider);
		g.setColor(Color.BLACK);
		g.drawString(title(), 20, height - 20);
		if (overlay && o != null)
		{
			System.out.println("here0");
			o.draw(g, width, height - 50, slider);
		}
	}
	
	abstract public void drawValues(Graphics g, int width, int height, double slider);
	
	abstract public String title();
	
	public Plot addOverlay(final Overlay newOverlay)
	{
		o = newOverlay;
		return this;
	}
}
