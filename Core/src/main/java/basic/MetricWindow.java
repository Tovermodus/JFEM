package basic;

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;
import java.awt.image.BufferedImage;
import java.util.Objects;
import java.util.concurrent.ConcurrentSkipListMap;
import java.util.concurrent.Executors;

public class MetricWindow
	extends JFrame
	implements KeyListener, WindowListener, ComponentListener, MouseWheelListener, MetricWindowInterface
{
	private final Canvas canvas;
	private final ConcurrentSkipListMap<String, Metric> metrics;
	private final DrawThread d;
	private volatile String currentPlot = null;
	private static MetricWindow INSTANCE;
	
	public static MetricWindowInterface getInstance()
	{
		try
		{
			if (GraphicsEnvironment.isHeadless())
				return new MetricWindowInterface()
				{
				};
		} catch (final HeadlessException e)
		{
			return new MetricWindowInterface()
			{
			};
		}
		if (GraphicsEnvironment.isHeadless())
			return new MetricWindowInterface()
			{
			};
		
		if (INSTANCE == null)
		{
			INSTANCE = new MetricWindow();
		}
		return INSTANCE;
	}
	
	private MetricWindow()
	{
		canvas = new Canvas();
		d = new DrawThread(canvas, canvas.getWidth(), canvas.getHeight());
		metrics = new ConcurrentSkipListMap<>();
		try
		{
			setSize(800, 800);
			setLayout(new BorderLayout());
			canvas.setFocusable(false);
			final JPanel pan = new JPanel();
			add(canvas, BorderLayout.CENTER);
			pan.setLayout(new BorderLayout());
			add(pan, BorderLayout.NORTH);
			pan.setFocusable(true);
			this.setFocusable(true);
			setVisible(true);
			addComponentListener(this);
			addKeyListener(this);
			canvas.addKeyListener(this);
			addMouseWheelListener(this);
			addWindowListener(this);
			Executors.newSingleThreadExecutor()
			         .execute(d);
		} catch (final HeadlessException e)
		{
		
		}
	}
	
	@Override
	public void addMetric(final Metric plot)
	{
		metrics.put("" + metrics.size(), plot);
		if (metrics.size() == 1)
		{
			d.setMetric(plot);
			currentPlot = metrics.firstKey();
		}
	}
	
	@Override
	public <M extends Metric> M setMetric(final String name, final M plot)
	{
		metrics.put(name, plot);
		if (metrics.size() == 1)
		{
			currentPlot = metrics.firstKey();
		}
		if (Objects.equals(currentPlot, name))
			d.setMetric(plot);
		return plot;
	}
	
	@Override
	public void componentResized(final ComponentEvent e)
	{
		d.setSize(canvas.getWidth(), canvas.getHeight());
	}
	
	@Override
	public void componentMoved(final ComponentEvent e)
	{
	}
	
	@Override
	public void componentShown(final ComponentEvent e)
	{
	}
	
	@Override
	public void componentHidden(final ComponentEvent e)
	{
	}
	
	@Override
	public void keyTyped(final KeyEvent e)
	{
	}
	
	@Override
	public void keyPressed(final KeyEvent e)
	{
		if (e.getKeyCode() == KeyEvent.VK_DOWN)
		{
			currentPlot = metrics.higherKey(currentPlot);
			if (currentPlot == null)
				currentPlot = metrics.firstKey();
		}
		if (e.getKeyCode() == KeyEvent.VK_UP)
		{
			currentPlot = metrics.lowerKey(currentPlot);
			if (currentPlot == null)
				currentPlot = metrics.lastKey();
		}
		if (metrics.size() != 0) d.setMetric(metrics.get(currentPlot));
	}
	
	@Override
	public void keyReleased(final KeyEvent e)
	{
	}
	
	@Override
	public void mouseWheelMoved(final MouseWheelEvent e)
	{
	}
	
	@Override
	public void windowOpened(final WindowEvent e)
	{
	}
	
	@Override
	public void windowClosing(final WindowEvent e)
	{
		d.running = false;
		synchronized (d)
		{
			try
			{
				System.out.println("requested stop, waiting");
				d.wait();
				System.out.println("drawthread stopped");
			} catch (final InterruptedException interruptedException)
			{
				interruptedException.printStackTrace();
			}
		}
		dispose();
	}
	
	@Override
	public void windowClosed(final WindowEvent e)
	{
		System.exit(0);
	}
	
	@Override
	public void windowIconified(final WindowEvent e)
	{
	}
	
	@Override
	public void windowDeiconified(final WindowEvent e)
	{
	}
	
	@Override
	public void windowActivated(final WindowEvent e)
	{
	}
	
	@Override
	public void windowDeactivated(final WindowEvent e)
	{
	}
	
	static class DrawThread
		implements Runnable
	{
		private volatile Metric currentPlot;
		volatile boolean running = true;
		volatile int width;
		volatile int height;
		volatile boolean isChecked;
		final Canvas canvas;
		BufferedImage content;
		
		public DrawThread(final Canvas canvas, final int width, final int height)
		{
			this.canvas = canvas;
			this.width = width;
			this.height = height;
			isChecked = false;
			content = new BufferedImage(Math.max(canvas.getWidth(), 1), Math.max(canvas.getHeight(), 1),
			                            BufferedImage.TYPE_INT_RGB);
		}
		
		@Override
		public void run()
		{
			while (running)
			{
				draw();
				try
				{
					Thread.sleep(30);
				} catch (final InterruptedException e)
				{
					e.printStackTrace();
				}
			}
			synchronized (this)
			{
				notifyAll();
				System.out.println("notified");
			}
			System.out.println("drawThread finished");
		}
		
		public synchronized void draw()
		{
			try
			{
				
				content.getGraphics()
				       .setColor(Color.white);
				content.getGraphics()
				       .fillRect(0, 0, width + 100, height + 100);
				if (currentPlot != null)
					currentPlot.draw(content.getGraphics(),
					                 width,
					                 height,
					                 isChecked);
				final Graphics g = content.getGraphics();
				g.setColor(Color.GREEN);
				String s = ((MetricWindow) getInstance()).currentPlot;
				if (s == null)
					s = "";
				g.drawString(s, width - 100, 100);
				canvas.getGraphics()
				      .drawImage(content, 0, 0, null);
			} catch (final Exception e)
			{
				e.printStackTrace();
			}
		}
		
		public synchronized void setMetric(final Metric p)
		{
			currentPlot = p;
		}
		
		public synchronized void setSize(final int width, final int height)
		{
			this.width = width;
			this.height = height;
			if (canvas.getWidth() > 0 && canvas.getHeight() > 0)
				content = new BufferedImage(canvas.getWidth(),
				                            canvas.getHeight(),
				                            BufferedImage.TYPE_INT_RGB);
		}
	}
}
