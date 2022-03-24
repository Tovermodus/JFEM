package dlm;

import io.vavr.Tuple2;
import linalg.DenseMatrix;
import linalg.DenseVector;

import javax.swing.*;
import java.awt.*;
import java.io.*;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.Executors;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public abstract class DLMSystemFrame
	extends JFrame
	implements DLMSystemInterface
{
	private final double dt;
	private final int timeSteps;
	final Fluid backGround;
	final List<Particle> particles;
	final DLMSolver solver;
	DenseMatrix fluidHistory;
	List<DenseMatrix> particleHistory;
	List<DenseMatrix> lagrangeHistory;
	volatile boolean running;
	volatile boolean shouldBeRunning;
	AtomicInteger iteration;
	final JButton store;
	JTextField it;
	final JButton load;
	final JButton showButton;
	final JButton stoprun;
	final String name;
	
	public DLMSystemFrame(final double dt,
	                      final int timeSteps,
	                      final Fluid backGround,
	                      final List<Particle> particles,
	                      final DLMSolver solver, final String name)
	{
		super("DLMActions");
		this.backGround = backGround;
		this.particles = particles;
		this.dt = dt;
		this.timeSteps = timeSteps;
		this.solver = solver;
		this.name = name;
		running = false;
		shouldBeRunning = false;
		iteration = new AtomicInteger(0);
		final GridLayout l = new GridLayout(3, 2);
		setLayout(l);
		stoprun = new JButton("currently: running. Click To stop");
		load = new JButton("load");
		store = new JButton("store");
		showButton = new JButton("show");
		it = new JTextField("Click to start");
		this.add(stoprun);
		this.add(load);
		this.add(store);
		this.add(it);
		this.add(showButton);
		this.setSize(400, 400);
		this.setVisible(true);
		setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		setRunText();
		stoprun.addActionListener(actionEvent ->
		                          {
			                          shouldBeRunning = !shouldBeRunning;
			                          setRunText();
			                          if (shouldBeRunning)
			                          {
				                          Executors.newSingleThreadExecutor()
				                                   .execute(() ->
				                                            {
					                                            if (iteration.get() == 0)
						                                            initLoop();
					                                            else
					                                            {
						                                            final var state
							                                            = readHistory(
							                                            iteration.get());
						                                            loop(iteration.get(),
						                                                 state._1,
						                                                 state._2);
					                                            }
				                                            });
			                          }
		                          });
		store.addActionListener(e ->
		                        {
			                        final File dir = new File("/home/tovermodus/dlm/" + name);
			                        if (!dir.mkdirs())
				                        throw new IllegalStateException("could not create dir");
			                        try
			                        {
				                        System.out.println(iteration.get());
				                        FileOutputStream fout = new FileOutputStream(
					                        "/home/tovermodus/dlm/" + name + "/fh");
				                        ObjectOutputStream outs = new ObjectOutputStream(fout);
				                        outs.writeObject(fluidHistory);
				                        outs.close();
				                        for (int i = 0; i < particles.size(); i++)
				                        {
					                        fout = new FileOutputStream(
						                        "/home/tovermodus/dlm/" + name + "/ph" + i);
					                        outs = new ObjectOutputStream(fout);
					                        outs.writeObject(particleHistory.get(i));
					                        outs.close();
					                        fout = new FileOutputStream(
						                        "/home/tovermodus/dlm/" + name + "/lh" + i);
					                        outs = new ObjectOutputStream(fout);
					                        outs.writeObject(lagrangeHistory.get(i));
					                        outs.close();
				                        }
				                        Runtime.getRuntime()
				                               .exec("cp -r /home/tovermodus/IdeaProjects/JFEM" +
					                                     "/Examples/src/dlm " +
					                                     "/home/tovermodus/dlm/" + name + "/code" +
					                                     "" + System.getProperty(
					                               "sun.java.command"));
			                        } catch (final FileNotFoundException ex)
			                        {
				                        ex.printStackTrace();
			                        } catch (final IOException ioException)
			                        {
				                        ioException.printStackTrace();
			                        }
		                        });
		showButton.addActionListener(e ->
		                             {
			                             Executors.newSingleThreadExecutor()
			                                      .execute(() ->
			                                               {
				                                               if (iteration.get() >= 1)
				                                               {
					                                               final var state
						                                               = readHistory(
						                                               iteration.get());
					                                               show(state._1,
					                                                    state._2,
					                                                    iteration.get());
				                                               }
			                                               });
		                             });
		load.addActionListener(e ->
		                       {
			                       try
			                       {
				                       FileInputStream fin = new FileInputStream(
					                       "/home/tovermodus/dlm/" + name + "/fh");
				                       ObjectInputStream ins = new ObjectInputStream(fin);
				                       fluidHistory = (DenseMatrix) ins.readObject();
				                       for (int i = 1; i < fluidHistory.getRows(); i++)
					                       if (fluidHistory.getRow(i)
					                                       .euclidianNorm() == 0)
					                       {
						                       iteration = new AtomicInteger(i - 1);
						                       System.out.println(iteration);
						                       break;
					                       }
				                       ins.close();
				                       particleHistory = new ArrayList<>();
				                       lagrangeHistory = new ArrayList<>();
				                       for (int i = 0; i < particles.size(); i++)
				                       {
					                       fin = new FileInputStream(
						                       "/home/tovermodus/dlm/" + name + "/ph" + i);
					                       ins = new ObjectInputStream(fin);
					                       particleHistory.add((DenseMatrix) ins.readObject());
					                       ins.close();
					                       fin = new FileInputStream(
						                       "/home/tovermodus/dlm/" + name + "/lh" + i);
					                       ins = new ObjectInputStream(fin);
					                       lagrangeHistory.add((DenseMatrix) ins.readObject());
					                       ins.close();
				                       }
			                       } catch (final FileNotFoundException ex)
			                       {
				                       ex.printStackTrace();
			                       } catch (final IOException ioException)
			                       {
				                       ioException.printStackTrace();
			                       } catch (final ClassNotFoundException classNotFoundException)
			                       {
				                       classNotFoundException.printStackTrace();
			                       }
		                       });
	}
	
	private void setRunText()
	{
		it.setText("Currently in iteration " + iteration.get() + " out of " + timeSteps);
		String text = "Currently: ";
		if (running)
		{
			showButton.setEnabled(false);
			//store.setEnabled(false);
			load.setEnabled(false);
			text += "running ";
		} else
		{
			showButton.setEnabled(true);
			//store.setEnabled(true);
			load.setEnabled(true);
			text += "stopped ";
		}
		if (running != shouldBeRunning)
		{
			stoprun.setEnabled(false);
			text += ". Waiting to be ";
			if (shouldBeRunning)
				text += "started";
			else
				text += "stopped";
		} else
		{
			stoprun.setEnabled(true);
		}
		stoprun.setText(text);
	}
	
	private void initLoop()
	{
		final double time = 0;
		final FluidIterate fluidState = backGround.buildInitialIterate();
		final List<ParticleIterate> particlestates = particles.stream()
		                                                      .map(p -> p.buildInitialIterate(dt))
		                                                      .collect(Collectors.toList());
		fluidHistory = new DenseMatrix(timeSteps + 1, fluidState.current.getLength());
		particleHistory = particlestates.stream()
		                                .map(p -> new DenseMatrix(timeSteps + 1, p.current.getLength()))
		                                .collect(Collectors.toList());
		lagrangeHistory = particlestates.stream()
		                                .map(p -> new DenseMatrix(timeSteps + 1, p.currentLagrange.getLength()))
		                                .collect(Collectors.toList());
		writeHistory(fluidState, particlestates, 0);
		postIterationCallback(fluidState, particlestates, time);
		iteration.set(0);
		loop(0, fluidState, particlestates);
	}
	
	private void loop(final int startStep, FluidIterate fluidState, List<ParticleIterate> particlestates)
	{
		running = true;
		setRunText();
		double time = startStep * dt;
		for (int i = startStep; i < timeSteps && shouldBeRunning; i++)
		{
			time += dt;
			final Tuple2<FluidIterate, List<ParticleIterate>> state = timeStep(fluidState,
			                                                                   particlestates,
			                                                                   time);
			fluidState = state._1;
			particlestates = state._2;
			postIterationCallback(fluidState, particlestates, time);
			writeHistory(fluidState, particlestates, i + 1);
			iteration.set(i + 1);
			setRunText();
		}
		running = false;
		setRunText();
	}
	
	private void writeHistory(final FluidIterate fluidState,
	                          final List<ParticleIterate> particlestates,
	                          final int step)
	{
		fluidHistory.addRow(fluidState.current, step);
		IntStream.range(0, particles.size())
		         .forEach(i -> particleHistory.get(i)
		                                      .addRow(particlestates.get(i).current, step));
		IntStream.range(0, particles.size())
		         .forEach(i -> lagrangeHistory.get(i)
		                                      .addRow(particlestates.get(i).currentLagrange, step));
	}
	
	private Tuple2<FluidIterate, List<ParticleIterate>> readHistory(final int step)
	{
		final DenseVector currentFluid = fluidHistory.getRow(step);
		final FluidIterate fluidIterate = new FluidIterate(currentFluid);
		final List<ParticleIterate> particleIterates = new ArrayList<>();
		for (int i = 0; i < particles.size(); i++)
		{
			final DenseVector currentParticle = particleHistory.get(i)
			                                                   .getRow(step);
			final DenseVector currentParticleLagrange = lagrangeHistory.get(i)
			                                                           .getRow(step);
			final DenseVector lastParticle = particleHistory.get(i)
			                                                .getRow(step - 1);
			final DenseVector lastParticleLagrange = lagrangeHistory.get(i)
			                                                        .getRow(step - 1);
			final ParticleIterate iterate = new ParticleIterate(currentParticle, currentParticleLagrange,
			                                                    lastParticle, lastParticleLagrange);
			particleIterates.add(iterate);
		}
		return new Tuple2<>(fluidIterate, particleIterates);
	}
	
	@Override
	public List<Particle> getParticles()
	{
		return particles;
	}
	
	@Override
	public Fluid getFluid()
	{
		return backGround;
	}
	
	@Override
	public DLMSolver getSolver()
	{
		return solver;
	}
	
	@Override
	public double getDT()
	{
		return dt;
	}
	
	@Override
	abstract public void postIterationCallback(final FluidIterate fluidState,
	                                           final List<ParticleIterate> particleStates,
	                                           final double time);
	
	@Override
	abstract public void show(final FluidIterate fluidState, final List<ParticleIterate> particleStates,
	                          final int iteration);
}
