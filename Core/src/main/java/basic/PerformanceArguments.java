package basic;

import java.util.Objects;

public class PerformanceArguments
{
	public static int GMRESResidual = 0;
	public static int GMRESGamma = 1;
	private static PerformanceArguments INSTANCE;
	public final boolean parallelizeThreads;
	public final int threadNumber;
	public final boolean executeChecks;
	public final int cacheSize;
	public final double doubleTolerance;
	public final int GMResData;
	
	private PerformanceArguments(final Boolean parallelizeThreads,
	                             final Integer threadNumber,
	                             final Boolean executeChecks,
	                             final Integer cacheSize,
	                             final Double doubleTolerance,
	                             final Integer GMResData)
	{
		this.parallelizeThreads = Objects.requireNonNullElse(parallelizeThreads, true);
		this.threadNumber = Objects.requireNonNullElse(threadNumber, 12);
		this.executeChecks = Objects.requireNonNullElse(executeChecks, true);
		this.cacheSize = Objects.requireNonNullElse(cacheSize, 50000);
		this.doubleTolerance = Objects.requireNonNullElse(doubleTolerance, 1e-10);
		this.GMResData = Objects.requireNonNullElse(GMResData, GMRESGamma);
	}
	
	public static void createInstance(final Boolean parallelizeThreads,
	                                  final Integer threadNumber,
	                                  final Boolean executeChecks,
	                                  final Integer cacheSize,
	                                  final Double doubleTolerance,
	                                  final Integer GMResData)
	{
		if (INSTANCE != null)
			throw new IllegalStateException();
		INSTANCE = new PerformanceArguments(parallelizeThreads, threadNumber, executeChecks, cacheSize,
		                                    doubleTolerance, GMResData);
	}
	
	public static PerformanceArguments getInstance()
	{
		if (INSTANCE == null)
			createInstance(null, null, null, null, null, null);
		return INSTANCE;
	}
	
	public static class PerformanceArgumentBuilder
	{
		public Boolean parallelizeThreads = null;
		public Integer threadNumber = null;
		public Boolean executeChecks = null;
		public Integer cacheSize = null;
		public Double doubleTolerance = null;
		public Integer GMResData = null;
		
		public void build()
		{
			PerformanceArguments.createInstance(parallelizeThreads, threadNumber, executeChecks,
			                                    cacheSize, doubleTolerance, GMResData);
		}
	}
}
