package basic;

import java.awt.*;
import java.util.ArrayList;

public abstract class Plot
{
	ArrayList<Overlay> overlays;
	
	public Plot()
	{
		overlays = new ArrayList<>();
	}
	
	void draw(final Graphics g, final int width, final int height, final double slider, final boolean overlay)
	{
		g.setColor(Color.BLACK);
		drawValues(g, width, height - 50, slider);
		g.setColor(Color.BLACK);
		g.drawString(title(), 20, height - 20);
		if (overlay)
		{
			for (final Overlay o : overlays)
				o.draw(g, width, height - 50, slider);
		}
	}
	
	abstract public void drawValues(Graphics g, int width, int height, double slider);
	
	abstract public String title();
	
	public Plot addOverlay(final Overlay newOverlay)
	{
		overlays.add(newOverlay);
		return this;
	}
}
