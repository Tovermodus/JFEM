package dlm;

import basic.PerformanceArguments;
import basic.PlotWindow;
import basic.StatLogger;
import distorted.DistortedOverlay;
import linalg.CoordinateVector;
import linalg.DenseMatrix;
import linalg.IntCoordinates;
import mixed.*;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentSkipListMap;

public class DLMFSIBenchmark
	extends DLMSystemFrame
{
	List<CoordinateVector> plotPoints;
	private Map<CoordinateVector, CoordinateVector> velocityValues;
	private Map<CoordinateVector, Double> pressureValues;
	
	public DLMFSIBenchmark(final DLMBenchmark.DLMBenchmarkConfig config,
	                       final int timeSteps,
	                       final MultiGridFluid backGround,
	                       final List<Particle> particles, final String name)
	{
		super(config.dt, timeSteps, backGround, particles, new DLMHybridMGSolver(config.smootherSteps,
		                                                                         config.overlap,
		                                                                         config.stepsCoarser,
		                                                                         backGround,
		                                                                         particles,
		                                                                         1.),
		      name);
		plotPoints = backGround.getSpace()
		                       .generatePlotPoints(41);
		velocityValues = new ConcurrentSkipListMap<>();
		pressureValues = new ConcurrentSkipListMap<>();
		System.out.println("Simulating up to time " + timeSteps * config.dt);
	}
	
	static double dt = 0.004;
	
	public static void main(final String[] args)
	{
		StatLogger.clear();
		final DLMBenchmark.DLMBenchmarkConfig config = new DLMBenchmark.DLMBenchmarkConfig();
		System.out.println(config);
		System.out.println(Runtime.getRuntime()
		                          .maxMemory());
		System.out.println(Runtime.getRuntime()
		                          .maxMemory() / 1000 + "KB");
		System.out.println(Runtime.getRuntime()
		                          .maxMemory() / 1000000 + "MB");
		System.out.println(Runtime.getRuntime()
		                          .maxMemory() / 1000000000 + "GB");
		StatLogger.log(config.toString());
		StatLogger.log("MAx MEm " + Runtime.getRuntime()
		                                   .maxMemory());
		final var builder = new PerformanceArguments.PerformanceArgumentBuilder();
		builder.build();
		final MultiGridFluid fluid = new BenchmarkFluid(CoordinateVector.fromValues(0, 0),
		                                                CoordinateVector.fromValues(2.2, 0.41),
		                                                new IntCoordinates(16, 4),
		                                                1,
		                                                config,
		                                                1,
		                                                1e-3);
		final List<Particle> particles = new ArrayList<>();
		particles.add(new FixedSolidParticle(CoordinateVector.fromValues(0.2, 0.2),
		                                     0.05,
		                                     3,
		                                     1,
		                                     config.lamb,
		                                     config.mu,
		                                     150));
		particles.add(new SolidBrickParticle(CoordinateVector.fromValues(0.4, 0.2),
		                                     0.4,
		                                     0.02,
		                                     1,
		                                     2,
		                                     CoordinateVector.fromValues(0, 0),
		                                     1000, 10, 10));
		final DLMFSIBenchmark system = new DLMFSIBenchmark(config,
		                                                   2000,
		                                                   fluid,
		                                                   particles,
		                                                   "benchfiner_" + fluid.refinements);
	}
	
	@Override
	public void postIterationCallback(final FluidIterate fluidState,
	                                  final List<ParticleIterate> particleStates,
	                                  final double time)
	{
		System.out.println("ITeration at time " + time + " is finished");
		final var points = backGround.getSpace()
		                             .generatePlotPoints(100);
		final var velos = new MixedTPFESpaceFunction<QkQkFunction>(backGround.getSpace()
		                                                                     .getShapeFunctionMap(),
		                                                           fluidState.current).getVelocityFunction()
		                                                                              .valuesInPoints(points);
		final var ps = new MixedTPFESpaceFunction<QkQkFunction>(backGround.getSpace()
		                                                                  .getShapeFunctionMap(),
		                                                        fluidState.current).getPressureFunction()
		                                                                           .valuesInPoints(points);
		
		final BufferedWriter writer;
		try
		{
			writer = new BufferedWriter(new FileWriter("../dlm/data " + String.format("%8.4e", time)));
			for (final var p : points)
			{
				final CoordinateVector vel = velos.get(p);
				final double pr = ps.get(p);
				final String str
					= String.format("%6.3e", p.x())
					+ "," + String.format("%6.3e", p.y())
					+ "," + String.format("%6.3e", vel.x())
					+ "," + String.format("%6.3e", vel.y())
					+ "," + String.format("%6.3e", pr) + "\n";
				writer.write(str);
			}
			writer.close();
		} catch (final IOException e)
		{
			e.printStackTrace();
		}
		final VelocityMagnitudePlot2D p =
			new VelocityMagnitudePlot2D(new MixedTPFESpaceFunction<>(backGround.getSpace()
			                                                                   .getShapeFunctionMap(),
			                                                         fluidState.current),
			                            backGround.getSpace()
			                                      .generatePlotPoints(41),
			                            41,
			                            "velocity and pressure" + time);
		for (int i = 0; i < particles.size(); i++)
		{
			final DenseMatrix particleHist = new DenseMatrix(1,
			                                                 particles.get(i)
			                                                          .getSystemSize());
			particleHist.addRow(particleStates.get(i).current, 0);
			final DistortedOverlay d = new DistortedOverlay(p,
			                                                particles.get(i)
			                                                         .getSpace(),
			                                                particleHist,
			                                                5);
			p.addOverlay(d);
		}
		PlotWindow.addPlotShow(p);
	}
	
	@Override
	public void show(final FluidIterate fluidState,
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
