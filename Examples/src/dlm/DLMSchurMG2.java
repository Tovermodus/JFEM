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
import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentSkipListMap;

public class DLMSchurMG2
	extends DLMSystem
{
	List<CoordinateVector> plotPoints;
	private final Map<CoordinateVector, CoordinateVector> velocityValues;
	private final Map<CoordinateVector, Double> pressureValues;
	
	public DLMSchurMG2(final double dt,
	                   final int timeSteps,
	                   final MultiGridFluid backGround,
	                   final List<Particle> particles)
	{
		super(dt, timeSteps, backGround, particles, new DLMFluidMGSolver(backGround));
		plotPoints = backGround.getSpace()
		                       .generatePlotPoints(30);
		velocityValues = new ConcurrentSkipListMap<>();
		pressureValues = new ConcurrentSkipListMap<>();
		System.out.println("Simulating up to time " + timeSteps * dt);
	}
	
	static double dt = 0.002;
	
	public static void main(final String[] args)
	{
		final var builder = new PerformanceArguments.PerformanceArgumentBuilder();
		builder.GMResData = PerformanceArguments.GMRESResidual;
		builder.build();
		final MultiGridFluid fluid = new BackgroundFluidMG(CoordinateVector.fromValues(0, 0),
		                                                   CoordinateVector.fromValues(1, 1),
		                                                   new IntCoordinates(4, 4),
		                                                   1,
		                                                   1,
		                                                   dt,
		                                                   0.5,
		                                                   5);
		final List<Particle> particles = new ArrayList<>();
		particles.add(new SolidParticle(CoordinateVector.fromValues(0.4, 0.5),
		                                0.05,
		                                2,
		                                1,
		                                CoordinateVector.fromValues(10, 3),
		                                1000,
		                                1000,
		                                15));
		particles.add(new SolidParticle(CoordinateVector.fromValues(0.6, 0.5),
		                                0.05,
		                                2,
		                                1,
		                                CoordinateVector.fromValues(-10, 3),
		                                1000,
		                                1000,
		                                15));
		particles.add(new Membrane(CoordinateVector.fromValues(0.5, 0.5),
		                           0.39,
		                           0.4,
		                           0,
		                           1,
		                           new CoordinateVector(2),
		                           10,
		                           10,
		                           1));
		final DLMSchurMG2 system = new DLMSchurMG2(dt, 10, fluid, particles);
		system.loop();
		system.summarize();
	}
	
	private void summarize()
	{
		final MixedPlot2DTime p = new MixedPlot2DTime(pressureValues,
		                                              velocityValues,
		                                              30,
		                                              "velocity and pressure");
		for (int i = 0; i < particles.size(); i++)
			p.addOverlay(particles.get(i)
			                      .generateOverlay(particleHistory.get(i), p));
		PlotWindow.addPlot(p);
	}
	
	@Override
	protected void postIterationCallback(final FluidIterate fluidState,
	                                     final List<ParticleIterate> particleStates,
	                                     final double time)
	{
		final MixedFunctionOnCells<TPCell, TPFace> velocityPressure
			= backGround.getVelocityPressure(fluidState);
		velocityValues.putAll(velocityPressure.velocityValuesInPointsAtTime(plotPoints, time));
		velocityValues.putAll(velocityPressure.velocityValuesInPointsAtTime(plotPoints, time));
		pressureValues.putAll(velocityPressure.pressureValuesInPointsAtTime(plotPoints, time));
		System.out.println("ITeration at time " + time + " is finished");
	}
}
