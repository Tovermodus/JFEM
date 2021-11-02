package basic;

import java.awt.*;

public interface Metric
{
	void draw(Graphics graphics, int width, int height, boolean isChecked);
	
	Metric getCopy();
}
