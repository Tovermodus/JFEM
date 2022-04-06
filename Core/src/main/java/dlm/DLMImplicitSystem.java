package dlm;

import basic.PlotWindow;
import distorted.DistortedOverlay;
import io.vavr.Tuple2;
import linalg.DenseMatrix;
import linalg.DenseVector;
import linalg.Vector;
import mixed.MixedTPFESpaceFunction;
import mixed.VelocityMagnitudePlot2D;
import org.jetbrains.annotations.NotNull;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;

public abstract class DLMImplicitSystem
	extends DLMSystemFrame
{
	public DLMImplicitSystem(final double dt,
	                         final int timeSteps,
	                         final Fluid backGround,
	                         final List<Particle> particles, final DLMSolver solver, final String name)
	{
		super(dt, timeSteps, backGround, particles, solver, name);
	}
	
	@Override
	public Tuple2<FluidIterate, List<ParticleIterate>> timeStep(final FluidIterate fluidState,
	                                                            final List<ParticleIterate> particleStates,
	                                                            final double time)
	{
		FluidIterate fluidGuess = new FluidIterate(fluidState.current);
		List<ParticleIterate> particleGuesses
			= particleStates.stream()
			                .map(p -> new ParticleIterate(p.current,
			                                              p.currentLagrange,
			                                              p.last,
			                                              p.lastLagrange))
			                .collect(Collectors.toList());
		for (int nonlinearIteration = 0; nonlinearIteration < 30; nonlinearIteration++)
		{
			final var system
				= createSystem(time, fluidState, particleStates, fluidGuess, particleGuesses);
			final DenseVector guessVector = new DenseVector(getTotalSystemSize());
			guessVector.addSmallVectorAt(fluidGuess.current, 0);
			int offset = fluidGuess.current.getLength();
			for (int j = 0; j < particles.size(); j++)
			{
				guessVector.addSmallVectorAt(particleGuesses.get(j).current, offset);
				offset += particles.get(j)
				                   .getSystemSize();
				guessVector.addSmallVectorAt(particleGuesses.get(j).currentLagrange, offset);
				offset += particles.get(j)
				                   .getLagrangeSize();
			}
			final DenseVector defect = new DenseVector(system._4.sub(system._3.mvMul(guessVector)));
			System.out.println("NONLIN DEFECT AFTER ITER" + nonlinearIteration + " TIME " + time + " NORM" +
				                   " " + defect.euclidianNorm() + " " + system._4.euclidianNorm());
			final Vector correction = getSolver().solve(system._3,
			                                            defect,
			                                            fluidGuess,
			                                            particleGuesses,
			                                            system._1,
			                                            system._2,
			                                            getDT(),
			                                            time);
			System.out.println("NONLIN Correct AFTER ITER" + nonlinearIteration + " TIME " + time + " " +
				                   "NORM" +
				                   " " + correction.euclidianNorm() + " " + guessVector.euclidianNorm());
			
			final Vector newGuess = guessVector.add(correction.mul(0.5));
			final DenseVector defect2 = new DenseVector(system._4.sub(system._3.mvMul(newGuess)));
			System.out.println("NONLIN DEFECT AFTER ITER" + nonlinearIteration + " TIME " + time + " NORM" +
				                   " " + defect2.euclidianNorm() + " " + system._4.euclidianNorm());
			final var newGuesses = createNewIterates(particleStates, newGuess);
			System.out.println(fluidGuess.current.sub(newGuesses._1.current)
			                                     .euclidianNorm() + " DIFFF");
			if (fluidGuess.current.sub(newGuesses._1.current)
			                      .euclidianNorm() < 1e-1)
				nonlinearIteration = 40;
			fluidGuess = newGuesses._1;
			particleGuesses = newGuesses._2;
			System.out.println();
			final VelocityMagnitudePlot2D p =
				new VelocityMagnitudePlot2D(
					new MixedTPFESpaceFunction<>(getFluid().getSpace()
					                                       .getShapeFunctionMap(),
					                             fluidGuess.current),
					getFluid().getSpace()
					          .generatePlotPoints(41),
					41, "-----");
			for (int i = 0; i < particles.size(); i++)
			{
				final DenseMatrix particleHist = new DenseMatrix(1,
				                                                 particles.get(i)
				                                                          .getSystemSize());
				particleHist.addRow(particleGuesses.get(i).current, 0);
				final DistortedOverlay d = new DistortedOverlay(p,
				                                                particles.get(i)
				                                                         .getSpace(),
				                                                particleHist,
				                                                5);
				p.addOverlay(d);
			}
			PlotWindow.addPlotShow(p);
		}
		final VelocityMagnitudePlot2D p =
			new VelocityMagnitudePlot2D(
				new MixedTPFESpaceFunction<>(getFluid().getSpace()
				                                       .getShapeFunctionMap(),
				                             fluidGuess.current),
				getFluid().getSpace()
				          .generatePlotPoints(81),
				81);
		for (int i = 0; i < particles.size(); i++)
		{
			final DenseMatrix particleHist = new DenseMatrix(1,
			                                                 particles.get(i)
			                                                          .getSystemSize());
			particleHist.addRow(particleGuesses.get(i).current, 0);
			final DistortedOverlay d = new DistortedOverlay(p,
			                                                particles.get(i)
			                                                         .getSpace(),
			                                                particleHist,
			                                                5);
			p.addOverlay(d);
		}
		PlotWindow.addPlotShow(p);
		return new Tuple2<>(fluidGuess, particleGuesses);
	}
	
	@Override
	@NotNull
	public Tuple2<FluidIterate, List<ParticleIterate>> createNewIterates(final List<ParticleIterate> previousGuess,
	                                                                     final Vector newGuess)
	{
		int offset;
		final FluidIterate ret = new FluidIterate(newGuess.slice(0, getFluid().getSystemSize()));
		offset = getFluid().getSystemSize();
		final List<ParticleIterate> iterates = new ArrayList<>();
		for (int i = 0; i < getParticles().size(); i++)
		{
			final ParticleIterate it
				= new ParticleIterate(newGuess.slice(offset,
				                                     offset + getParticles().get(i)
				                                                            .getSystemSize())
				                              .mul(getDT()),
				                      newGuess.slice(offset + getParticles().get(i)
				                                                            .getSystemSize(),
				                                     offset + getParticles().get(i)
				                                                            .getSystemSize()
					                                     + getParticles().get(i)
					                                                     .getLagrangeSize()),
				                      previousGuess.get(i).last, previousGuess.get(i).lastLagrange);
			iterates.add(it);
			offset += getParticles().get(i)
			                        .getSystemSize() + getParticles().get(i)
			                                                         .getLagrangeSize();
		}
		return new Tuple2<>(ret, iterates);
	}
}
