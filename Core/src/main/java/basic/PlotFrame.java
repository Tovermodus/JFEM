package basic;

import basic.ScalarFunction;
import com.google.common.collect.Iterables;
import linalg.CoordinateVector;
import linalg.Vector;

import java.awt.*;
import java.awt.event.*;
import java.awt.image.BufferedImage;
import java.awt.image.ImageObserver;
import java.util.*;
import java.util.List;
import java.util.concurrent.*;

public class PlotFrame
{
	CoordinateVector startCoordinates;
	CoordinateVector endCoordinates;
	Map<String, Map<CoordinateVector, Double>> valueList;
	Iterator<String> current;
	
	String currentTitle;
	int units;
	int dimension;
	ScheduledExecutorService executorService;
	boolean drawing = false;
	private void drawCurrentValues(Graphics gr, int width, int height)
	{
		if(!drawing)
		{
			drawing = true;
			Map<CoordinateVector, Double> currentValues = valueList.get(currentTitle);
			int gridPointsPerSpaceDim = (int) (Math.pow(currentValues.size(), 1. / dimension));
			List<CoordinateVector> drawnVectors = new ArrayList<>(currentValues.keySet());
			double currentZ =
				startCoordinates.at(dimension - 1) + units / 100. * (endCoordinates.at(dimension - 1) - startCoordinates.at(dimension - 1));
			if (dimension == 3)
			{
				drawnVectors.sort(Comparator.comparingDouble(v -> -Math.abs(v.z() - currentZ)));
				System.out.println(drawnVectors.get(0).z() + " " + Iterables.getLast(drawnVectors).z() + " " + currentZ);
				drawnVectors = drawnVectors.subList(drawnVectors.size()-(int)(3*gridPointsPerSpaceDim*gridPointsPerSpaceDim),
					drawnVectors.size()-1);
			}
			
			double[] range = getRange();
			BufferedImage img = new BufferedImage(2000,2000,BufferedImage.TYPE_INT_RGB);
			Graphics g = img.createGraphics();
			g.setColor(Color.white);
			g.fillRect(0,0,2000,2000);
			drawnVectors.forEach(vector ->
			{
				
				// .getHeight()-100)/gridPointsPerSpaceDim+1);
				g.setColor(convertToColor(currentValues.get(vector),range));
				Vector relativeCoords = vector.sub(startCoordinates);
				int posx =
					50 + (int) (width * relativeCoords.at(0) / (endCoordinates.at(0) - startCoordinates.at(0)));
				int posy =
					50 + (int) (height * relativeCoords.at(1) / (endCoordinates.at(1) - startCoordinates.at(1)));
				g.fillRect(posx, posy, width / gridPointsPerSpaceDim + 2, height / gridPointsPerSpaceDim + 2);
			});
			g.setColor(Color.BLACK);
			g.drawString("title: " + currentTitle + " max: " + range[1] + " min: " + range[0], 40, 40);
			if(dimension ==3)
			{
				g.drawLine(40, 50, 40, 50 + height);
				g.drawLine(30, 50, 50, 50);
				g.drawLine(30, 50 + height, 50, 50 + height);
				g.drawLine(30, 50 + height, 50, 50 + height);
				g.drawLine(30,
					50 + (int) (height *  units/100.), 50, 50 + (int) (height * units/100.));
			}
			
			gr.drawImage(img, 0, 0, null);
			drawing = false;
		}
	}
	private double[] getRange()
	{
		Map<CoordinateVector, Double> currentValues = valueList.get(currentTitle);
		OptionalDouble maxopt =
			currentValues.values().stream().mapToDouble(Double::doubleValue).max();
		double max = 0;
		if(maxopt.isPresent())
			max = maxopt.getAsDouble();
		OptionalDouble minopt =
			currentValues.values().stream().mapToDouble(Double::doubleValue).min();
		double min = 0;
		if(minopt.isPresent())
			min = minopt.getAsDouble();
		return new double[]{min, max};
	}
	private Color convertToColor(double value,double[] range)
	{
		double v =
			( value - range[0])/(range[1] - range[0]);
		int c = (int)(255*v);
		int gr = 0;
		if(c>=255)
		{
			gr = c - 255;
			c = 255;
		}
		if(c <= 0)
		{
			gr = -c;
			c = 0;
		}
		return new Color(c,gr,255-c);
	}
	public PlotFrame(List<Map<CoordinateVector, Double>> valueMapList, CoordinateVector startCoordinates,
	                 CoordinateVector endCoordinates)
	{
		Map<String,Map<CoordinateVector, Double>> vmap= new TreeMap<>();
		for(int i = 0; i < valueMapList.size(); i++)
			vmap.put("Values "+i, valueMapList.get(i));
		initialize(vmap,startCoordinates,endCoordinates);
	}
	public PlotFrame(List<ScalarFunction> functions, List<CoordinateVector> points,
	                 CoordinateVector startCoordinates,
	                 CoordinateVector endCoordinates)
	{
		Map<String,Map<CoordinateVector, Double>> vmap= new TreeMap<>();
		for(int i = 0; i < functions.size(); i++)
			vmap.put("Values "+i, functions.get(i).valuesInPoints(points));
		initialize(vmap,startCoordinates,endCoordinates);
	}
	public PlotFrame(Map<String,Map<CoordinateVector, Double>> valueList, CoordinateVector startCoordinates,
	                 CoordinateVector endCoordinates)
	{
		initialize(valueList, startCoordinates, endCoordinates);
	}
	public PlotFrame(Map<String,Map<CoordinateVector, Double>> valueList, CoordinateVector startCoordinates,
	                 CoordinateVector endCoordinates, double startTime, double endTime)
	{
		initialize(valueList, startCoordinates, endCoordinates);
	}
	private void initialize(Map<String,Map<CoordinateVector, Double>> valueList, CoordinateVector startCoordinates,
	                        CoordinateVector endCoordinates)
	{
		this.startCoordinates = startCoordinates;
		this.endCoordinates = endCoordinates;
		this.valueList = valueList;
		executorService = Executors.newSingleThreadScheduledExecutor();
		
		dimension = startCoordinates.getLength();
		current = this.valueList.keySet().iterator();
		currentTitle = current.next();
		units = 0;
		Frame j = new Frame("plot");
		int size = 1000;
		double maxx = endCoordinates.at(0);
		double maxy = endCoordinates.at(1);
		double minx = startCoordinates.at(0);
		double miny = startCoordinates.at(1);
		j.setBounds(2100,200,
			(int)(size*(maxx-minx)),(int)(size*(maxy-miny)));
		j.setVisible(true);
		j.addKeyListener(new KeyAdapter()
		{
			@Override
			public void keyPressed(KeyEvent e)
			{
				if(e.getKeyCode() == KeyEvent.VK_UP)
					units--;
				if(e.getKeyCode() == KeyEvent.VK_DOWN)
					units++;
				if(e.getKeyCode() == KeyEvent.VK_RIGHT)
				{
					if(current.hasNext())
						currentTitle = current.next();
					else
						current = valueList.keySet().iterator();
				}
				System.out.println("KKKKKK");
				if (units < 0)
					units += 100;
				if (units > 100)
					units -= 100;
				System.out.println(units);
				Graphics g = j.getGraphics();
				//g.clearRect(0,0,2000,2000);
				executorService.schedule(new Runnable()
				{
					@Override
					public void run()
					{
						drawCurrentValues(j.getGraphics(), j.getWidth() - 300,
							j.getHeight() - 300);
					}
				},0,TimeUnit.SECONDS);
				
				j.setVisible(true);
			}
		});
		j.addMouseWheelListener(new MouseWheelListener()
		{
			@Override
			public void mouseWheelMoved(MouseWheelEvent e)
			{
				units += e.getWheelRotation();
				if (units < 0)
					units += 100;
				if (units > 100)
					units -= 100;
				System.out.println(units);
				Graphics g = j.getGraphics();
				//g.clearRect(0,0,2000,2000);
					executorService.schedule(new Runnable()
					{
						@Override
						public void run()
						{
							drawCurrentValues(j.getGraphics(), j.getWidth() - 300,
								j.getHeight() - 300);
						}
					},0,TimeUnit.SECONDS);
				
				j.setVisible(true);
			}
		});
		j.addMouseListener(new MouseAdapter()
		                   {
			                   int i = 0;
			                   @Override
			                   public void mouseClicked(MouseEvent e)
			                   {
				                   Graphics g = j.getGraphics();
				                   g.clearRect(0,0,2000,2000);
				                   if(current.hasNext())
				                   	currentTitle = current.next();
				                   else
				                   	current = valueList.keySet().iterator();
				                   executorService.schedule(new Runnable()
				                   {
					                   @Override
					                   public void run()
					                   {
						                   drawCurrentValues(j.getGraphics(), j.getWidth() - 300,
							                   j.getHeight() - 300);
					                   }
				                   },0,  TimeUnit.SECONDS);
				                   if(i >= valueList.size())
				                   	i = 0;
				                   j.setVisible(true);
			                   }
			
			
			                   @Override
			                   public void mouseExited(MouseEvent e)
			                   {
				                   j.dispose();
				                   executorService.shutdown();
				                   System.exit(0);
			                   }
		                   }
		
		);
	}
}
