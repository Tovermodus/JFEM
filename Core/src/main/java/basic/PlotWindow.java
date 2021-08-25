package basic;

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;
import java.awt.image.BufferedImage;
import java.util.concurrent.CopyOnWriteArrayList;
import java.util.concurrent.Executors;

public class PlotWindow extends JFrame implements KeyListener, WindowListener, ComponentListener, MouseWheelListener
{
	final JSlider slider;
	final JCheckBox overlay;
	final Canvas canvas;
	final CopyOnWriteArrayList<Plot> plots;
	final DrawThread d;
	int currentPlot = 0;
	
	public PlotWindow()
	{
		
		setSize(800, 800);
		setLayout(new BorderLayout());
		slider = new JSlider(JSlider.HORIZONTAL, 0, 100, 0);
		slider.setFocusable(false);
		overlay = new JCheckBox("Show Overlay");
		canvas = new Canvas();
		canvas.setFocusable(false);
		final JPanel pan = new JPanel();
		add(canvas, BorderLayout.CENTER);
		pan.setLayout(new BorderLayout());
		pan.add(slider, BorderLayout.CENTER);
		pan.add(overlay, BorderLayout.EAST);
		add(pan, BorderLayout.NORTH);
		pan.setFocusable(true);
		this.setFocusable(true);
		setVisible(true);
		d = new DrawThread(canvas, canvas.getWidth(), canvas.getHeight(), slider);
		plots = new CopyOnWriteArrayList<>();
		addComponentListener(this);
		addKeyListener(this);
		canvas.addKeyListener(this);
		addMouseWheelListener(this);
		addWindowListener(this);
		overlay.addChangeListener(e ->
		                          {
			                          d.isChecked = overlay.isSelected();
		                          });
		Executors.newSingleThreadExecutor().execute(d);
	}
	
	public void addPlot(final Plot plot)
	{
		if (plots.size() == 0) d.setPlot(plot);
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
		int newValue = slider.getValue();
		if (e.getKeyCode() == KeyEvent.VK_DOWN) currentPlot--;
		if (e.getKeyCode() == KeyEvent.VK_UP) currentPlot++;
		if (e.getKeyCode() == KeyEvent.VK_RIGHT) newValue++;
		if (e.getKeyCode() == KeyEvent.VK_LEFT) newValue--;
		if (newValue < 0) newValue = 99;
		if (newValue > 99) newValue = 0;
		if (plots.size() != 0) currentPlot = currentPlot % plots.size();
		while (currentPlot < 0) currentPlot += plots.size();
		System.out.println(currentPlot);
		if (plots.size() != 0) d.setPlot(plots.get(currentPlot));
		slider.setValue(newValue);
	}
	
	@Override
	public void keyReleased(final KeyEvent e)
	{
	}
	
	@Override
	public void mouseWheelMoved(final MouseWheelEvent e)
	{
		int newValue = (slider.getValue() + e.getWheelRotation());
		if (newValue < 0) newValue = 99;
		slider.setValue(newValue % 100);
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
}

class DrawThread implements Runnable
{
	private volatile Plot currentPlot;
	volatile boolean running = true;
	volatile int width;
	volatile int height;
	private final JSlider slider;
	volatile boolean isChecked;
	final Canvas canvas;
	BufferedImage content;
	
	public DrawThread(final Canvas canvas, final int width, final int height, final JSlider slider)
	{
		this.canvas = canvas;
		this.width = width;
		this.height = height;
		this.slider = slider;
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
		content.getGraphics().setColor(Color.white);
		content.getGraphics().fillRect(0, 0, width + 100, height + 100);
		if (currentPlot != null)
			currentPlot.draw(content.getGraphics(), width, height, 0.01 * slider.getValue(), isChecked);
		canvas.getGraphics().drawImage(content, 0, 0, null);
	}
	
	public synchronized void setPlot(final Plot p)
	{
		currentPlot = p;
	}
	
	public synchronized void setSize(final int width, final int height)
	{
		this.width = width;
		this.height = height;
		content = new BufferedImage(canvas.getWidth(), canvas.getHeight(), BufferedImage.TYPE_INT_RGB);
	}
}
