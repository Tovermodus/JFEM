package basic;

public interface PlotWindowInterface
{
	default DrawThreadInterface getDrawThread()
	{
		return new DrawThreadInterface();
	}
	
	default void addPlotPrivate(final Plot plot)
	{
	}
	
	class DrawThreadInterface
	{
		public void setPlot(final Plot plot)
		{
		}
	}
}
