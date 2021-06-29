package basic;

import java.util.Objects;

public class PerformanceArguments
{
	private static PerformanceArguments INSTANCE;
	public final boolean parallelizeThreads;
	public final int threadNumber;
	public final boolean executeChecks;
	public final int cacheSize;
	
	private PerformanceArguments(Boolean parallelizeThreads, Integer threadNumber, Boolean executeChecks,
	                             Integer cacheSize) {
		this.parallelizeThreads = Objects.requireNonNullElse(parallelizeThreads, true);
		this.threadNumber = Objects.requireNonNullElse(threadNumber, 12);
		this.executeChecks = Objects.requireNonNullElse(executeChecks, false);
		this.cacheSize = Objects.requireNonNullElse(cacheSize, 50000);
	}
	public static void createInstance(Boolean parallelizeThreads, Integer threadNumber, Boolean executeChecks,
	                                  Integer cacheSize) {
		INSTANCE = new PerformanceArguments(parallelizeThreads, threadNumber, executeChecks, cacheSize);
	}
	public static PerformanceArguments getInstance() {
		if(INSTANCE == null)
			createInstance(null, null, null, null);
		return INSTANCE;
	}
	
	public static class PerformanceArgumentBuilder{
		public Boolean parallelizeThreads = null;
		public Integer threadNumber = null;
		public Boolean executeChecks = null;
		public Integer cacheSize = null;
		public void build()
		{
			PerformanceArguments.createInstance(parallelizeThreads, threadNumber, executeChecks, cacheSize);
		}
	}
}
