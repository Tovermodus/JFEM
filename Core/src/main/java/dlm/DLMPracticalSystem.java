package dlm;

import io.vavr.Tuple2;
import it.unimi.dsi.fastutil.ints.Int2DoubleArrayMap;
import it.unimi.dsi.fastutil.ints.Int2DoubleMap;
import linalg.*;

import javax.swing.*;
import java.awt.*;
import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.concurrent.Executors;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public abstract class DLMPracticalSystem
	extends JFrame
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
	
	public DLMPracticalSystem(final double dt,
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
	
	private Tuple2<FluidIterate, List<ParticleIterate>> timeStep(final FluidIterate fluidState,
	                                                             final List<ParticleIterate> particleStates,
	                                                             final double time)
	{
		final Map<IntCoordinates, SparseMatrix> blocks = new HashMap<>();
		final DenseVector rhs =
			new DenseVector(backGround.getSystemSize() + particles.stream()
			                                                      .mapToInt(p -> p.getSystemSize() + p.getLagrangeSize())
			                                                      .sum());
		final FluidSystem fluidSystem = backGround.buildSystem(time, fluidState, particles);
		final List<ParticleSystem> particleSystems =
			IntStream.range(0, particles.size())
			         .mapToObj(i -> particles.get(i)
			                                 .buildSystem(backGround, time, particleStates.get(i)))
			         .collect(Collectors.toList());
		
		int offset = 0;
		final var fluidBlockRhs = Fluid.getBlockRhs(fluidSystem, dt);
		blocks.put(new IntCoordinates(0, 0), fluidBlockRhs._1);
		rhs.addSmallVectorAt(fluidBlockRhs._2, 0);
		offset += fluidBlockRhs._1.getCols();
		for (int i = 0; i < particles.size(); i++)
		{
			offset = addParticleBlocks(blocks, rhs, particleSystems, offset, i);
		}
		final BlockSparseMatrix systemMatrix = new BlockSparseMatrix(blocks, rhs.getLength(), rhs.getLength());
		final Tuple2<BlockSparseMatrix, DenseVector> system = applyBoundaryValues(systemMatrix, rhs, time);
		final Vector solution = solver.solve(system._1, system._2, fluidState, particleStates, fluidSystem,
		                                     particleSystems, dt, time);
		final FluidIterate ret = new FluidIterate(solution.slice(0, backGround.getSystemSize()));
		offset = backGround.getSystemSize();
		final List<ParticleIterate> iterates = new ArrayList<>();
		for (int i = 0; i < particles.size(); i++)
		{
			final ParticleIterate it
				= new ParticleIterate(particleStates.get(i),
				                      solution.slice(offset,
				                                     offset + particles.get(i)
				                                                       .getSystemSize())
				                              .mul(dt),
				                      solution.slice(offset + particles.get(i)
				                                                       .getSystemSize(),
				                                     offset + particles.get(i)
				                                                       .getSystemSize()
					                                     + particles.get(i)
					                                                .getLagrangeSize()));
			iterates.add(it);
			offset += particles.get(i)
			                   .getSystemSize() + particles.get(i)
			                                               .getLagrangeSize();
		}
		return new Tuple2<>(ret, iterates);
	}
	
	protected Tuple2<BlockSparseMatrix, DenseVector> applyBoundaryValues(final BlockSparseMatrix systemMatrix,
	                                                                     final DenseVector rhs, final double t)
	{
		final SparseMatrix ret = systemMatrix.toSparse();
		final DenseVector retRhs = new DenseVector(rhs);
		final Int2DoubleMap nodeValues = new Int2DoubleArrayMap();
		final Int2DoubleMap fluidDirichletNodeValues = backGround.getDirichletNodeValues(t);
		nodeValues.putAll(fluidDirichletNodeValues);
		for (int i = 0; i < particles.size(); i++)
		{
			final Int2DoubleMap particleDirichletNodeValues = particles.get(i)
			                                                           .getDirichletNodeValues(t);
			final int finalI = i;
			final int particleSpaceSize = particles.get(i)
			                                       .getSystemSize();
			particleDirichletNodeValues.forEach((node, val) ->
				                                    nodeValues.put(node + systemMatrix.getBlockStarts()[finalI + 1],
				                                                   val.doubleValue()));
			particleDirichletNodeValues.forEach((node, val) ->
				                                    nodeValues.put(node + particleSpaceSize + systemMatrix.getBlockStarts()[finalI + 1],
				                                                   val.doubleValue()));
		}
		nodeValues.forEach((node, val) ->
		                   {
			                   final DenseVector column = ret.getColumn(node);
			                   for (int i = 0; i < column.getLength(); i++)
			                   {
				                   retRhs.add(-column.at(i) * val, i);
			                   }
			                   ret.deleteColumn(node);
			                   ret.deleteRow(node);
			                   ret.set(1, node, node);
			                   retRhs.set(val, node);
		                   });
		return new Tuple2<>(new BlockSparseMatrix(ret, systemMatrix.getBlockStarts()), retRhs);
	}
	
	private int addParticleBlocks(final Map<IntCoordinates, SparseMatrix> blocks,
	                              final DenseVector rhs,
	                              final List<ParticleSystem> particleSystems,
	                              int offset,
	                              final int i)
	{
		final var particleBlockRhs = particles.get(i)
		                                      .getBlockRhs(particleSystems.get(i), dt);
		blocks.put(new IntCoordinates(offset, offset), particleBlockRhs._1);
		rhs.addSmallVectorAt(particleBlockRhs._2, offset);
		final var particleBackgroundLagrangeBlock =
			particles.get(i)
			         .getLagrangeBackgroundBlock(particleSystems.get(i), backGround, dt);
		blocks.put(new IntCoordinates(0, offset), particleBackgroundLagrangeBlock);
		final var particleBackgroundLagrangeBlockTranspose =
			particles.get(i)
			         .getLagrangeBackgroundBlockTranspose(particleSystems.get(i), backGround, dt);
		blocks.put(new IntCoordinates(offset, 0), particleBackgroundLagrangeBlockTranspose);
		offset += particleBlockRhs._1.getCols();
		return offset;
	}
	
	protected abstract void postIterationCallback(final FluidIterate fluidState,
	                                              final List<ParticleIterate> particleStates,
	                                              final double time);
	
	protected abstract void show(final FluidIterate fluidState,
	                             final List<ParticleIterate> particleStates,
	                             int iteration);
}
