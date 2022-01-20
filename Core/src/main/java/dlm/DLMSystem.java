package dlm;

import io.vavr.Tuple2;
import linalg.*;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public abstract class DLMSystem
{
	private final double dt;
	private final int timeSteps;
	Fluid backGround;
	List<Particle> particles;
	DenseMatrix fluidHistory;
	List<DenseMatrix> particleHistory;
	List<DenseMatrix> lagrangeHistory;
	
	public DLMSystem(final double dt, final int timeSteps, final Fluid backGround, final List<Particle> particles)
	{
		this.backGround = backGround;
		this.particles = particles;
		this.dt = dt;
		this.timeSteps = timeSteps;
	}
	
	public void loop()
	{
		double time = 0;
		FluidIterate fluidState = backGround.buildInitialIterate();
		List<ParticleIterate> particlestates = particles.stream()
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
		for (int i = 0; i < timeSteps; i++)
		{
			time += dt;
			final Tuple2<FluidIterate, List<ParticleIterate>> state = timeStep(fluidState,
			                                                                   particlestates,
			                                                                   time);
			fluidState = state._1;
			particlestates = state._2;
			postIterationCallback(fluidState, particlestates, time);
			writeHistory(fluidState, particlestates, i + 1);
		}
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
	
	private Tuple2<FluidIterate, List<ParticleIterate>> timeStep(final FluidIterate fluidState,
	                                                             final List<ParticleIterate> particleStates,
	                                                             final double time)
	{
		final Map<IntCoordinates, SparseMatrix> blocks = new HashMap<>();
		final DenseVector rhs =
			new DenseVector(backGround.getSystemSize() + particles.stream()
			                                                      .mapToInt(p -> p.getSystemSize() + p.getLagrangeSize())
			                                                      .sum());
		final FluidSystem fluidSystem = backGround.buildSystem(time, fluidState);
		final List<ParticleSystem> particleSystems =
			IntStream.range(0, particles.size())
			         .mapToObj(i -> particles.get(i)
			                                 .buildSystem(backGround, time, particleStates.get(i)))
			         .collect(Collectors.toList());
		
		int offset = 0;
		final var fluidBlockRhs = backGround.getBlockRhs(fluidSystem, dt);
		blocks.put(new IntCoordinates(0, 0), fluidBlockRhs._1);
		rhs.addSmallVectorAt(fluidBlockRhs._2, 0);
		offset += fluidBlockRhs._1.getCols();
		for (int i = 0; i < particles.size(); i++)
		{
			offset = addParticleBlocks(blocks, rhs, particleSystems, offset, i);
		}
		final BlockSparseMatrix systemMatrix = new BlockSparseMatrix(blocks, rhs.getLength(), rhs.getLength());
		final DenseVector solution = solve(systemMatrix, rhs);
		final FluidIterate ret = new FluidIterate(solution.slice(0, backGround.getSystemSize()));
		offset = backGround.getSystemSize();
		final List<ParticleIterate> iterates = new ArrayList<>();
		for (int i = 0; i < particles.size(); i++)
		{
			final ParticleIterate it
				= new ParticleIterate(particleStates.get(i),
				                      solution.slice(offset,
				                                     offset + particles.get(i)
				                                                       .getSystemSize()),
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
	
	protected abstract DenseVector solve(BlockSparseMatrix systemMatrix, DenseVector rhs);
	
	protected abstract void postIterationCallback(final FluidIterate fluidState,
	                                              final List<ParticleIterate> particleStates,
	                                              final double time);
}
