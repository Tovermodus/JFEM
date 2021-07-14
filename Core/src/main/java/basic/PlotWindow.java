package basic;

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;
import java.awt.image.BufferedImage;
import java.util.ArrayList;
import java.util.concurrent.*;

public class PlotWindow extends JFrame implements KeyListener, WindowListener, ComponentListener, MouseWheelListener
{
	final JSlider slider;
	final Canvas canvas;
	final CopyOnWriteArrayList<Plot> plots;
	final DrawThread d;
	int currentPlot = 0;
	public PlotWindow() {
		
		setSize(800,800);
		setLayout(new BorderLayout());
		slider = new JSlider(JSlider.HORIZONTAL, 0, 100,0);
		slider.setFocusable(false);
		canvas = new Canvas();
		canvas.setFocusable(false);
		add(slider, BorderLayout.NORTH);
		add(canvas, BorderLayout.CENTER);
		setVisible(true);
		plots = new CopyOnWriteArrayList<>();
		addComponentListener(this);
		addKeyListener(this);
		addMouseWheelListener(this);
		addWindowListener(this);
		d = new DrawThread(canvas, canvas.getWidth(), canvas.getHeight(), slider);
		Executors.newSingleThreadExecutor().execute(d);
		
	}
	public void addPlot(Plot plot)
	{
		if(plots.size() == 0)
			d.setPlot(plot);
		plots.add(plot);
	}
	
	@Override
	public void componentResized(ComponentEvent e)
	{
		d.setSize(canvas.getWidth(), canvas.getHeight());
	}
	@Override
	public void componentMoved(ComponentEvent e)
	{
	}
	
	@Override
	public void componentShown(ComponentEvent e)
	{
	}
	
	@Override
	public void componentHidden(ComponentEvent e)
	{
	}
	
	@Override
	public void keyTyped(KeyEvent e)
	{
	}
	
	@Override
	public void keyPressed(KeyEvent e)
	{
		int newValue = slider.getValue();
		if(e.getKeyCode() == KeyEvent.VK_DOWN)
			currentPlot--;
		if(e.getKeyCode() == KeyEvent.VK_UP)
			currentPlot++;
		if(e.getKeyCode() == KeyEvent.VK_RIGHT)
			newValue++;
		if(e.getKeyCode() == KeyEvent.VK_LEFT)
			newValue--;
		if(newValue < 0)
			newValue = 99;
		if(newValue > 99)
			newValue = 0;
		currentPlot = currentPlot % plots.size();
		System.out.println(currentPlot);
		d.setPlot(plots.get(currentPlot));
		slider.setValue(newValue);
	}
	
	@Override
	public void keyReleased(KeyEvent e)
	{
	}
	
	@Override
	public void mouseWheelMoved(MouseWheelEvent e)
	{
		int newValue = (slider.getValue() + e.getWheelRotation());
		if(newValue < 0)
			newValue = 99;
		slider.setValue(newValue%100);
	}
	
	@Override
	public void windowOpened(WindowEvent e)
	{
	}
	
	@Override
	public void windowClosing(WindowEvent e)
	{
		d.running = false;
		synchronized (d)
		{
			try
			{
				System.out.println("requested stop, waiting");
				d.wait();
				System.out.println("drawthread stopped");
			} catch (InterruptedException interruptedException)
			{
				interruptedException.printStackTrace();
			}
		}
		dispose();
	}
	
	@Override
	public void windowClosed(WindowEvent e)
	{
		System.exit(0);
	}
	
	@Override
	public void windowIconified(WindowEvent e)
	{
	}
	
	@Override
	public void windowDeiconified(WindowEvent e)
	{
	}
	
	@Override
	public void windowActivated(WindowEvent e)
	{
	}
	
	@Override
	public void windowDeactivated(WindowEvent e)
	{
	}
}

class DrawThread implements Runnable {
	private volatile Plot currentPlot;
	volatile boolean running = true;
	volatile int width;
	volatile int height;
	private final JSlider slider;
	final Canvas canvas;
	BufferedImage content;
	
	public DrawThread(Canvas canvas, int width, int height, JSlider slider)
	{
		this.canvas = canvas;
		this.width = width;
		this.height = height;
		this.slider = slider;
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
			} catch (InterruptedException e)
			{
				e.printStackTrace();
			}
		}
		synchronized(this)
		{
			notifyAll();
			System.out.println("notified");
		}
		System.out.println("drawThread finished");
	}
	public synchronized void draw()
	{
		content.getGraphics().setColor(Color.white);
		content.getGraphics().fillRect(0,0,width+100,height+100);
		if(currentPlot != null)
			currentPlot.drawValues(content.getGraphics(), width, height, 0.01*slider.getValue());
		canvas.getGraphics().drawImage(content,0,0,null);
	}
	public synchronized void setPlot(Plot p)
	{
		currentPlot = p;
	}
	public synchronized void setSize(int width, int height)
	{
		this.width = width;
		this.height = height;
		content = new BufferedImage(canvas.getWidth(), canvas.getHeight(), BufferedImage.TYPE_INT_RGB);
	}
}
