package dlm;

import basic.PlotWindow;
import distorted.DistortedVectorFESpaceFunction;
import linalg.BlockSparseMatrix;
import linalg.CoordinateVector;
import linalg.DenseVector;
import linalg.IntCoordinates;
import mixed.MixedFunctionOnCells;
import mixed.MixedPlot2DTime;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentSkipListMap;

public class DLM2ParticleSystem
	extends DLMSystem
{
	List<CoordinateVector> plotPoints;
	private final Map<CoordinateVector, CoordinateVector> velocityValues;
	private final Map<CoordinateVector, CoordinateVector> XPrimeValues;
	private final Map<CoordinateVector, Double> pressureValues;
	
	public DLM2ParticleSystem(final double dt,
	                          final int timeSteps,
	                          final Fluid backGround,
	                          final List<Particle> particles)
	{
		super(dt, timeSteps, backGround, particles);
		plotPoints = backGround.getSpace()
		                       .generatePlotPoints(30);
		velocityValues = new ConcurrentSkipListMap<>();
		XPrimeValues = new ConcurrentSkipListMap<>();
		pressureValues = new ConcurrentSkipListMap<>();
		System.out.println("Simulating up to time " + timeSteps * dt);
	}
	
	public static void main(final String[] args)
	{
		final Fluid fluid = new BackgroundFluid(CoordinateVector.fromValues(0, 0),
		                                        CoordinateVector.fromValues(1, 1),
		                                        new IntCoordinates(8, 8),
		                                        1,
		                                        1,
		                                        10);
		final List<Particle> particles = new ArrayList<>();
		particles.add(new SolidParticle(CoordinateVector.fromValues(0.2, 0.5),
		                                0.05,
		                                1,
		                                1,
		                                CoordinateVector.getUnitVector(2, 0)
		                                                .mul(0.0),
		                                1000,
		                                1000,
		                                10));
//		particles.add(new SolidParticle(CoordinateVector.fromValues(0.6, 0.5),
//		                                0.05,
//		                                1,
//		                                1,
//		                                CoordinateVector.getUnitVector(2, 0)
//		                                                .mul(-0.001),
//		                                1000,
//		                                1000,
//		                                1));
//		particles.add(new Membrane(CoordinateVector.fromValues(0.5, 0.5),
//		                           0.3,
//		                           0.4,
//		                           0,
//		                           1,
//		                           new CoordinateVector(2),
//		                           1000,
//		                           1000,
//		                           1));
		final DLM2ParticleSystem system = new DLM2ParticleSystem(0.001, 30, fluid, particles);
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
		System.out.println(particleHistory.get(0));
		System.out.println(lagrangeHistory.get(0));
	}
	
	@Override
	protected DenseVector solve(final BlockSparseMatrix systemMatrix, final DenseVector rhs)
	{
		return systemMatrix.toSparse()
		                   .solve(rhs);
	}
	
	@Override
	protected void postIterationCallback(final FluidIterate fluidState,
	                                     final List<ParticleIterate> particleStates,
	                                     final double time)
	{
		final MixedFunctionOnCells<TPCell, TPFace> velocityPressure
			= backGround.getVelocityPressure(fluidState);
		velocityValues.putAll(velocityPressure.velocityValuesInPointsAtTime(plotPoints, time));
		final DistortedVectorFESpaceFunction XPrime =
			new DistortedVectorFESpaceFunction(particles.get(0)
			                                            .getSpace()
			                                            .getShapeFunctions(),
			                                   particleStates.get(0).current.sub(particleStates.get(0).last)
			                                                                .mul(1. / 0.01));
		
		XPrimeValues.putAll(XPrime.valuesInPointsAtTime(plotPoints, time));
		velocityValues.putAll(velocityPressure.velocityValuesInPointsAtTime(plotPoints, time));
		pressureValues.putAll(velocityPressure.pressureValuesInPointsAtTime(plotPoints, time));
		System.out.println("ITeration at time " + time + " is finished");
	}
}
