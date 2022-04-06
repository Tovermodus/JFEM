package basic;

public interface MetricWindowInterface
{
	default void addMetric(final Metric plot)
	{
	}
	
	default <M extends Metric> M setMetric(final String name, final M plot)
	{
		return plot;
	}
}
