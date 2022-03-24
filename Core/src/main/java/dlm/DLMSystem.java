package dlm;

import io.vavr.Tuple2;
import it.unimi.dsi.fastutil.ints.Int2DoubleArrayMap;
import it.unimi.dsi.fastutil.ints.Int2DoubleMap;
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
	final Fluid backGround;
	final List<Particle> particles;
	final DLMSolver solver;
	DenseMatrix fluidHistory;
	List<DenseMatrix> particleHistory;
	List<DenseMatrix> lagrangeHistory;
	
	public DLMSystem(final double dt,
	                 final int timeSteps,
	                 final Fluid backGround,
	                 final List<Particle> particles,
	                 final DLMSolver solver)
	{
		this.backGround = backGround;
		this.particles = particles;
		this.dt = dt;
		this.timeSteps = timeSteps;
		this.solver = solver;
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
		final FluidSystem fluidSystem = backGround.buildSystem(time, fluidState, fluidState, particles);
		final List<ParticleSystem> particleSystems =
			IntStream.range(0, particles.size())
			         .mapToObj(i -> particles.get(i)
			                                 .buildSystem(backGround, time, particleStates.get(i),
			                                              particleStates.get(i)))
			         .collect(Collectors.toList());
		
		int offset = 0;
		final var fluidBlockRhs = Fluid.assembleBlockRhs(fluidSystem, dt);
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
			System.out.println("scaleparticle " + it.current.absMaxElement());
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
}
