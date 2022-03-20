package dlm;

import basic.PerformanceArguments;
import basic.PlotWindow;
import linalg.CoordinateVector;
import linalg.IntCoordinates;
import mixed.MixedFunctionOnCells;
import mixed.MixedPlot2DTime;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentSkipListMap;

public class DLMBenchmark
	extends DLMPracticalSystem
{
	List<CoordinateVector> plotPoints;
	private Map<CoordinateVector, CoordinateVector> velocityValues;
	private Map<CoordinateVector, Double> pressureValues;
	
	public DLMBenchmark(final double dt,
	                    final int timeSteps,
	                    final MultiGridFluid backGround,
	                    final List<Particle> particles, final String name)
	{
		super(dt, timeSteps, backGround, particles, new DLMHybridMGSolver(3, 3, backGround, particles), name);
		plotPoints = backGround.getSpace()
		                       .generatePlotPoints(41);
		velocityValues = new ConcurrentSkipListMap<>();
		pressureValues = new ConcurrentSkipListMap<>();
		System.out.println("Simulating up to time " + timeSteps * dt);
	}
	
	static double dt = 0.001;
	
	public static void main(final String[] args)
	{
		final var builder = new PerformanceArguments.PerformanceArgumentBuilder();
		builder.build();
		final MultiGridFluid fluid = new BenchmarkFluid(CoordinateVector.fromValues(0, 0),
		                                                CoordinateVector.fromValues(1.5, 0.41),
		                                                new IntCoordinates(16, 4),
		                                                1,
		                                                2,
		                                                dt,
		                                                1,
		                                                1e-3);
		final List<Particle> particles = new ArrayList<>();
		particles.add(new FixedSolidParticle(CoordinateVector.fromValues(0.15, 0.2),
		                                     0.05,
		                                     2,
		                                     1,
		                                     10000,
		                                     10000,
		                                     150));
		final DLMBenchmark system = new DLMBenchmark(dt,
		                                             2000,
		                                             fluid,
		                                             particles,
		                                             "benchhigherreyn" + fluid.refinements);
	}
	
	@Override
	protected void postIterationCallback(final FluidIterate fluidState,
	                                     final List<ParticleIterate> particleStates,
	                                     final double time)
	{
		System.out.println("ITeration at time " + time + " is finished");
	}
	
	@Override
	protected void show(final FluidIterate fluidState,
	                    final List<ParticleIterate> particleStates,
	                    final int iteration)
	{
		velocityValues = new HashMap<>();
		pressureValues = new HashMap<>();
		final int len = iteration + 1;
		final int[] indices;
		if (iteration > 100)
		{
			indices = new int[100];
			for (int i = 0; i < 100; i++)
				indices[i] = (int) (len / 99.0 * i);
		} else
		{
			indices = new int[iteration];
			for (int i = 0; i < iteration; i++)
				indices[i] = i;
		}
		for (int j = 0; j < indices.length; j++)
		{
			final int i = indices[j];
			final MixedFunctionOnCells<TPCell, TPFace> velocityPressure
				= backGround.getVelocityPressure(new FluidIterate(fluidHistory.getRow(i)));
			velocityValues.putAll(velocityPressure.velocityValuesInPointsAtTime(plotPoints, i * dt));
			pressureValues.putAll(velocityPressure.pressureValuesInPointsAtTime(plotPoints, i * dt));
			System.out.println("plot at time " + i * dt + " is finished");
		}
		final MixedPlot2DTime p = new MixedPlot2DTime(pressureValues,
		                                              velocityValues,
		                                              41,
		                                              "velocity and pressure");
		for (int i = 0; i < particles.size(); i++)
			p.addOverlay(
				particles.get(i)
				         .generateOverlay(
					         particleHistory.get(i)
					                        .slice(new IntCoordinates(0, 0),
					                               new IntCoordinates(iteration + 1,
					                                                  particles.get(i)
					                                                           .getSystemSize())),
					         p));
		PlotWindow.addPlot(p);
	}
}
