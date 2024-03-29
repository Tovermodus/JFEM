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

public class DLMRestartable
	extends DLMSystemFrame
{
	List<CoordinateVector> plotPoints;
	private Map<CoordinateVector, CoordinateVector> velocityValues;
	private Map<CoordinateVector, Double> pressureValues;
	
	public DLMRestartable(final double dt,
	                      final int timeSteps,
	                      final MultiGridFluid backGround,
	                      final List<Particle> particles, final String name)
	{
		super(dt,
		      timeSteps,
		      backGround,
		      particles,
		      //new DLMHybridMGSolver(3, 3, 3, backGround, particles, 1),
		      new DLMFluidMGSolver(backGround),
		      name);
		plotPoints = backGround.getSpace()
		                       .generatePlotPoints(41);
		velocityValues = new ConcurrentSkipListMap<>();
		pressureValues = new ConcurrentSkipListMap<>();
		System.out.println("Simulating up to time " + timeSteps * dt);
	}
	
	static double dt = 0.004;
	
	public static void main(final String[] args)
	{
		final var builder = new PerformanceArguments.PerformanceArgumentBuilder();
		builder.GMResData = PerformanceArguments.GMRESResidual;
		builder.build();
		final MultiGridFluid fluid
			= new BenchmarkFluid(CoordinateVector.fromValues(0, 0),
			                     CoordinateVector.fromValues(2, 1),
			                     new IntCoordinates(8, 4),
			                     1,
			                     new DLMBenchmark.DLMBenchmarkConfig(2,
			                                                         3,
			                                                         0.5,
			                                                         3,
			                                                         1,
			                                                         1,
			                                                         0.004,
			                                                         3),
			                     1,
			                     1e-3);
//			new BackgroundFluidMG(CoordinateVector.fromValues(0, 0),
//		                                                   CoordinateVector.fromValues(2, 1),
//		                                                   new IntCoordinates(8, 4),
//		                                                   1,
//		                                                   3,
//		                                                   dt,
//		                                                   1,
//		                                                   5);
		final List<Particle> particles = new ArrayList<>();
		particles.add(new SolidBrickParticle(CoordinateVector.fromValues(0.6, 0.45),
		                                     0.1,
		                                     0.45,
		                                     1,
		                                     2,
		                                     CoordinateVector.fromValues(0, 0),
		                                     10,
		                                     10,
		                                     5));
		particles.add(new SolidParticle(CoordinateVector.fromValues(0.75, 0.6),
		                                0.05,
		                                2,
		                                2,
		                                CoordinateVector.fromValues(0, 0),
		                                10,
		                                10,
		                                5));
		particles.add(new SolidParticle(CoordinateVector.fromValues(0.4, 0.5),
		                                0.05,
		                                2,
		                                2,
		                                CoordinateVector.fromValues(0, 0),
		                                10,
		                                10,
		                                5));
		particles.add(new Membrane(CoordinateVector.fromValues(0.6, 0.5),
		                           0.34,
		                           0.35,
		                           1,
		                           2,
		                           new CoordinateVector(2),
		                           10,
		                           10,
		                           5));
		particles.add(new FixedSolidParticle(CoordinateVector.fromValues(0.15, 0.25),
		                                     0.05,
		                                     2,
		                                     2,
		                                     100,
		                                     100,
		                                     5));
		particles.add(new FixedSolidParticle(CoordinateVector.fromValues(0.15, 0.75),
		                                     0.05,
		                                     2,
		                                     2,
		                                     100,
		                                     100,
		                                     5));
		particles.add(new FixedSolidParticle(CoordinateVector.fromValues(1.1, 0.25),
		                                     0.05,
		                                     2,
		                                     2,
		                                     100,
		                                     100,
		                                     5));
		particles.add(new FixedSolidParticle(CoordinateVector.fromValues(1.1, 0.75),
		                                     0.05,
		                                     2,
		                                     2,
		                                     100,
		                                     100,
		                                     5));
		final DLMRestartable system = new DLMRestartable(dt,
		                                                 400,
		                                                 fluid,
		                                                 particles,
		                                                 "tough" + fluid.refinements);
	}
	
	@Override
	public void postIterationCallback(final FluidIterate fluidState,
	                                  final List<ParticleIterate> particleStates,
	                                  final double time)
	{
		System.out.println("ITeration at time " + time + " is finished");
	}
	
	@Override
	public void show(final FluidIterate fluidState,
	                 final List<ParticleIterate> particleStates,
	                 final int iteration)
	{
		velocityValues = new HashMap<>();
		pressureValues = new HashMap<>();
		for (int i = 0; i < iteration + 1; i++)
		{
			
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
