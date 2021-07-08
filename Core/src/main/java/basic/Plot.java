package basic;

import javax.swing.*;
import java.awt.*;

public interface Plot
{
	default void draw(Graphics g, int width, int height, double slider)
	{
		g.setColor(Color.BLACK);
		drawValues(g,width,height-200,slider);
		g.drawString(title(), 20,height - 20);
	}
	void drawValues(Graphics g, int width, int height, double slider);
	String title();
}
