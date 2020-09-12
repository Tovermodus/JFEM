import linalg.CoordinateVector;
import linalg.Vector;

import java.awt.*;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.util.List;
import java.util.Map;
import java.util.OptionalDouble;

public class PlotFrame
{
	public PlotFrame(List<Map<CoordinateVector, Double>> valueList, CoordinateVector startCoordinates, CoordinateVector endCoordinates)
	{
		
		Frame j = new Frame("plot");
		int size = 1000;
		double maxx = endCoordinates.at(0);
		double maxy = endCoordinates.at(1);
		double minx = startCoordinates.at(0);
		double miny = startCoordinates.at(1);
		j.setBounds(2100,200,
			(int)(size*(maxx-minx)),(int)(size*(maxy-miny)));
		j.setVisible(true);
		j.addMouseListener(new MouseListener()
		                   {
			                   int i = 0;
			                   @Override
			                   public void mouseClicked(MouseEvent e)
			                   {
				                   Graphics g = j.getGraphics();
				                   g.clearRect(0,0,2000,2000);
				                   Map<CoordinateVector, Double> vals = valueList.get(i++);
				                   int n = (int)(Math.sqrt(vals.size()));
				                   OptionalDouble maxopt =
					                   vals.values().stream().mapToDouble(Double::doubleValue).max();
				                   double max = 0;
				                   if(maxopt.isPresent())
					                   max = maxopt.getAsDouble();
				                   OptionalDouble minopt = vals.values().stream().mapToDouble(Double::doubleValue).min();
				                   double min = 0;
				                   if(minopt.isPresent())
					                   min = minopt.getAsDouble();
				                   double finalMin = min;
				                   double finalMax = max;
				                   vals.keySet().forEach(coordinateVector -> {
					                   // .getHeight()-100)/n+1);
					                   double v;
					                   if(i == 1)
						                   v =
							                   ( vals.get(coordinateVector) - finalMin)/(finalMax - finalMin);
					                   else
						                   v = (vals.get(coordinateVector)+1.5)/3;
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
					                   g.setColor(new Color(c,gr,
						                   255-c));
					                   Vector relc = coordinateVector.sub(startCoordinates);
					                   int posx =
						                   50+(int)((j.getWidth()-100)*relc.at(0)/(endCoordinates.at(0) - startCoordinates.at(0)));
					                   int posy =
						                   50+(int)((j.getHeight()-100)*relc.at(1)/(endCoordinates.at(1) - startCoordinates.at(1)));
					                   g.fillRect(posx,posy,(j.getWidth()-100)/n+2,(j.getHeight()-100)/n+2);
				                   });
				                   g.setColor(Color.BLACK);
				                   g.drawString("i: "+i+" max: "+ max +" min: "+ min,40,40);
				                   if(i >= valueList.size())
				                   	i = 0;
				                   j.setVisible(true);
			                   }
			
			                   @Override
			                   public void mousePressed(MouseEvent e)
			                   {
				
			                   }
			
			                   @Override
			                   public void mouseReleased(MouseEvent e)
			                   {
				
			                   }
			
			                   @Override
			                   public void mouseEntered(MouseEvent e)
			                   {
				
			                   }
			
			                   @Override
			                   public void mouseExited(MouseEvent e)
			                   {
				                   System.exit(0);
			                   }
		                   }
		
		);
	}
}
