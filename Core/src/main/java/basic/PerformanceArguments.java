package basic;

import java.util.Objects;

public class PerformanceArguments
{
	private static PerformanceArguments INSTANCE;
	public final boolean parallelizeThreads;
	public final int threadNumber;
	public final boolean executeChecks;
	public final int cacheSize;
	public final double doubleTolerance;
	
	private PerformanceArguments(Boolean parallelizeThreads, Integer threadNumber, Boolean executeChecks,
	                             Integer cacheSize, Double doubleTolerance) {
		this.parallelizeThreads = Objects.requireNonNullElse(parallelizeThreads, true);
		this.threadNumber = Objects.requireNonNullElse(threadNumber, 12);
		this.executeChecks = Objects.requireNonNullElse(executeChecks, true);
		this.cacheSize = Objects.requireNonNullElse(cacheSize, 50000);
		this.doubleTolerance = Objects.requireNonNullElse(doubleTolerance, 1e-14);
	}
	public static void createInstance(Boolean parallelizeThreads, Integer threadNumber, Boolean executeChecks,
	                                  Integer cacheSize, Double doubleTolerance) {
		if(INSTANCE != null)
			throw new IllegalStateException();
		INSTANCE = new PerformanceArguments(parallelizeThreads, threadNumber, executeChecks, cacheSize, doubleTolerance);
	}
	public static PerformanceArguments getInstance() {
		if(INSTANCE == null)
			createInstance(null, null, null, null, null);
		return INSTANCE;
	}
	
	public static class PerformanceArgumentBuilder{
		public Boolean parallelizeThreads = null;
		public Integer threadNumber = null;
		public Boolean executeChecks = null;
		public Integer cacheSize = null;
		public Double doubleTolerance = null;
		public void build()
		{
			PerformanceArguments.createInstance(parallelizeThreads, threadNumber, executeChecks,
				cacheSize, doubleTolerance);
		}
	}
}
