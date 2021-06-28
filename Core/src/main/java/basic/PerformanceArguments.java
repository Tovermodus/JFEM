package basic;

import java.util.Objects;

public class PerformanceArguments
{
	private static PerformanceArguments INSTANCE;
	public final boolean parallelizeThreads;
	public final int threadNumber;
	public final boolean executeChecks;
	private PerformanceArguments(Boolean parallelizeThreads, Integer threadNumber, Boolean executeChecks) {
		this.parallelizeThreads = Objects.requireNonNullElse(parallelizeThreads, true);
		this.threadNumber = Objects.requireNonNullElse(threadNumber, 12);
		this.executeChecks = Objects.requireNonNullElse(executeChecks, false);
	}
	public static void createInstance(Boolean parallelizeThreads, Integer threadNumber, Boolean executeChecks) {
		INSTANCE = new PerformanceArguments(parallelizeThreads, threadNumber, executeChecks);
	}
	public static PerformanceArguments getInstance() {
		if(INSTANCE == null)
			createInstance(null, null, null);
		return INSTANCE;
	}
}
