package linalg;

import basic.Metric;

import java.awt.*;
import java.util.concurrent.CopyOnWriteArrayList;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.function.Supplier;

public class IterativeSolverConvergenceMetric
	implements Metric
{
	final int offset = 50;
	final int textOffset = 20;
	private CopyOnWriteArrayList<Double> residuals;
	public double goal;
	public int drawnPoints = 500;
	ExecutorService executorService;
	
	public IterativeSolverConvergenceMetric(final double goal)
	{
		this.goal = goal;
		residuals = new CopyOnWriteArrayList<>();
		executorService = Executors.newFixedThreadPool(2);
	}
	
	public synchronized void publishIterate(final double residual)
	{
		residuals.add(residual);
	}
	
	public void publishIterateAsync(final Supplier<Double> generatingFunction)
	{
		final IterativeSolverConvergenceMetric metric = this;
		final Runnable r = () ->
		{
			final Double residual = generatingFunction.get();
			metric.publishIterate(residual);
		};
		executorService.execute(r);
		//r.run();
	}
	
	public synchronized void restart()
	{
		if (residuals.size() > 0)
			residuals = new CopyOnWriteArrayList<>();
	}
	
	@Override
	public void draw(final Graphics graphics, final int width, final int height, final boolean isChecked)
	{
		final int size = (int) (Math.min(width, height) / drawnPoints + 3);
		
		final double maxResidual = residuals.stream()
		                                    .mapToDouble(x -> x)
		                                    .max()
		                                    .orElse(goal * 1000);
		graphics.setColor(Color.BLACK);
		graphics.fillRect(0, 0, width, height);
		graphics.setColor(Color.white);
		graphics.drawLine(offset, height - offset, offset, offset);
		graphics.drawLine(offset, height - offset, width - offset, height - offset);
		drawYMark(graphics, goal, height, maxResidual);
		drawYMark(graphics, maxResidual, height, maxResidual);
		drawYMark(graphics, 1, height, maxResidual);
		final int maxLog = (int) Math.log10(maxResidual);
		final int goalLog = (int) Math.log10(goal);
		for (int i = goalLog; i < maxLog; i++)
		{
			drawYMark(graphics, Math.pow(10, i), height, maxResidual);
		}
		final int xTicks = Math.max(1, (int) Math.pow(10, (int) Math.log10(residuals.size() * 1.2)) / 2);
		for (int i = 0; i <= residuals.size() / xTicks; i++)
		{
			drawXMark(graphics, i * xTicks, width, height);
		}
		drawXMark(graphics, residuals.size(), width, height);
		graphics.setColor(Color.BLUE);
		for (int i = 0; i < drawnPoints; i++)
		{
			final double x = offset + (width - 2.0 * offset) * i / drawnPoints;
			final double y = valueToY(pointToResidual(i), height, maxResidual);
			graphics.fillOval((int) x - size, (int) y - size, 2 * size, 2 * size);
			//graphics.drawLine((int) x + size, (int) y - size, (int) x - size, (int) y + size);
		}
	}
	
	@Override
	public Metric getCopy()
	{
		final IterativeSolverConvergenceMetric m = new IterativeSolverConvergenceMetric(goal);
		m.residuals = residuals;
		m.drawnPoints = drawnPoints;
		return m;
	}
	
	private void drawYMark(final Graphics g, final double value, final int height, final double max)
	{
		g.drawLine(offset - 5, (int) valueToY(value, height, max), offset + 5,
		           (int) valueToY(value, height, max));
		g.drawString("" + String.format("%2.1e", value), textOffset, (int) valueToY(value, height, max));
	}
	
	private void drawXMark(final Graphics g, final int iterate, final int width, final int height)
	{
		final int residSize = Math.max(1, residuals.size());
		final int x = (width - 2 * offset) * iterate / residSize;
		g.drawLine(offset + x,
		           height - offset - 5,
		           offset + x,
		           height - offset + 5);
		g.drawString("" + iterate,
		             offset + x - 10, height - textOffset);
	}
	
	private double valueToY(final double value, final int height, final double max)
	{
		return 1.0 * height - offset - (height * 1.0 - 2.0 * offset) * (Math.log(value) - Math.log(
			goal)) / (Math.log(
			max) - Math.log(goal));
	}
	
	private synchronized double pointToResidual(final int point)
	{
		final int nIterates = residuals.size();
		if (nIterates == 0)
			return goal;
		if (point == drawnPoints - 1)
			return residuals.get(nIterates - 1);
		final double val = 100;
		final int pointIndex = (int) (1.0 * nIterates * point / (drawnPoints + 1));
		final int pointIndex2 = (int) (1.0 * nIterates * (point + 1) / (drawnPoints + 1));
		return residuals.get(pointIndex);
	}
}
