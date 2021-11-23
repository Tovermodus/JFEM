package basic;

import java.awt.*;
import java.util.concurrent.atomic.AtomicInteger;

public class CountMetric
	implements Metric
{
	final int size;
	AtomicInteger i;
	
	public CountMetric(final int size)
	{
		this.size = size;
		i = new AtomicInteger(0);
	}
	
	@Override
	public void draw(final Graphics graphics, final int width, final int height, final boolean isChecked)
	{
		graphics.setColor(Color.BLACK);
		graphics.drawString("Iteration " + i.get() + " out of " + size + ", " + (int) (i.get() * 100 / size) + "%",
		                    50,
		                    50);
	}
	
	@Override
	public Metric getCopy()
	{
		throw new UnsupportedOperationException("askjdh");
	}
	
	public void reset()
	{
		i = new AtomicInteger(0);
	}
	
	public void increment()
	{
		i.incrementAndGet();
	}
}
