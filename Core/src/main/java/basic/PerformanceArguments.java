package basic;

import org.apache.log4j.Level;
import org.apache.log4j.Logger;
import org.apache.spark.SparkConf;
import org.apache.spark.api.java.JavaSparkContext;

import java.util.Objects;

public class PerformanceArguments
{
	private static PerformanceArguments INSTANCE;
	public final boolean parallelizeThreads;
	public final int threadNumber;
	public final boolean executeChecks;
	public final int cacheSize;
	public final JavaSparkContext sparkContext;
	public final SparkConf conf;
	
	private PerformanceArguments(Boolean parallelizeThreads, Integer threadNumber, Boolean executeChecks,
	                             Integer cacheSize) {
		this.parallelizeThreads = Objects.requireNonNullElse(parallelizeThreads, true);
		this.threadNumber = Objects.requireNonNullElse(threadNumber, 12);
		this.executeChecks = Objects.requireNonNullElse(executeChecks, true);
		this.cacheSize = Objects.requireNonNullElse(cacheSize, 50000);
		
		this.conf = new SparkConf().setAppName("JFEM").setMaster("local[*]");
		Logger.getLogger("org.apache.spark").setLevel(Level.WARN);
		sparkContext = new JavaSparkContext(conf);
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
