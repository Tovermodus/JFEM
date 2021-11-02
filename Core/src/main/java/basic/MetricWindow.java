package basic;

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;
import java.awt.image.BufferedImage;
import java.util.concurrent.CopyOnWriteArrayList;
import java.util.concurrent.Executors;

public class MetricWindow
	extends JFrame
	implements KeyListener, WindowListener, ComponentListener, MouseWheelListener
{
	private final Canvas canvas;
	private final CopyOnWriteArrayList<Metric> plots;
	private final DrawThread d;
	private int currentPlot = 0;
	private static MetricWindow INSTANCE;
	
	public static MetricWindow getInstance()
	{
		if (INSTANCE == null)
			INSTANCE = new MetricWindow();
		return INSTANCE;
	}
	
	private MetricWindow()
	{
		
		setSize(800, 800);
		setLayout(new BorderLayout());
		canvas = new Canvas();
		canvas.setFocusable(false);
		final JPanel pan = new JPanel();
		add(canvas, BorderLayout.CENTER);
		pan.setLayout(new BorderLayout());
		add(pan, BorderLayout.NORTH);
		pan.setFocusable(true);
		this.setFocusable(true);
		setVisible(true);
		d = new DrawThread(canvas, canvas.getWidth(), canvas.getHeight());
		plots = new CopyOnWriteArrayList<>();
		addComponentListener(this);
		addKeyListener(this);
		canvas.addKeyListener(this);
		addMouseWheelListener(this);
		addWindowListener(this);
		Executors.newSingleThreadExecutor()
		         .execute(d);
	}
	
	public void addMetric(final Metric plot)
	{
		if (plots.size() == 0) d.setMetric(plot);
		plots.add(plot);
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
		if (e.getKeyCode() == KeyEvent.VK_DOWN) currentPlot++;
		if (e.getKeyCode() == KeyEvent.VK_UP) currentPlot--;
		if (plots.size() != 0) currentPlot = currentPlot % plots.size();
		while (currentPlot < 0) currentPlot += plots.size();
		System.out.println(currentPlot);
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
			content = new BufferedImage(canvas.getWidth(), canvas.getHeight(), BufferedImage.TYPE_INT_RGB);
		}
		
		@Override
		public void run()
		{
			while (running)
			{
				draw();
				try
				{
					Thread.sleep(10);
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
			content.getGraphics()
			       .setColor(Color.white);
			content.getGraphics()
			       .fillRect(0, 0, width + 100, height + 100);
			if (currentPlot != null)
				currentPlot.draw(content.getGraphics(),
				                 width,
				                 height,
				                 isChecked);
			canvas.getGraphics()
			      .drawImage(content, 0, 0, null);
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
