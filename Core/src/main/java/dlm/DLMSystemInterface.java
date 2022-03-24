package dlm;

import io.vavr.Tuple2;
import io.vavr.Tuple4;
import it.unimi.dsi.fastutil.ints.Int2DoubleArrayMap;
import it.unimi.dsi.fastutil.ints.Int2DoubleMap;
import linalg.*;
import org.jetbrains.annotations.NotNull;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public interface DLMSystemInterface
{
	
	Fluid getFluid();
	
	List<Particle> getParticles();
	
	DLMSolver getSolver();
	
	double getDT();
	
	default Tuple2<FluidIterate, List<ParticleIterate>> timeStep(final FluidIterate fluidState,
	                                                             final List<ParticleIterate> particleStates,
	                                                             final double time)
	{
		
		final var system = createSystem(time, fluidState, particleStates, fluidState, particleStates);
		
		final Vector solution = getSolver().solve(system._3, system._4, fluidState, particleStates, system._1,
		                                          system._2, getDT(), time);
		
		return createNewIterates(particleStates, solution);
	}
	
	default Tuple4<FluidSystem, List<ParticleSystem>, BlockSparseMatrix, DenseVector>
	createSystem(final double time,
	             final FluidIterate fluidState,
	             final List<ParticleIterate> particleStates,
	             final FluidIterate fluidGuess,
	             final List<ParticleIterate> particleGuesses)
	{
		
		final Map<IntCoordinates, SparseMatrix> blocks = new HashMap<>();
		final DenseVector rhs = new DenseVector(getTotalSystemSize());
		final FluidSystem fluidSystem = getFluid().buildSystem(time, fluidState, fluidGuess, getParticles());
		final List<ParticleSystem> particleSystems =
			IntStream.range(0, getParticles().size())
			         .mapToObj(i -> getParticles().get(i)
			                                      .buildSystem(getFluid(),
			                                                   time,
			                                                   particleStates.get(i),
			                                                   particleGuesses.get(i)))
			         .collect(Collectors.toList());
		
		int offset = 0;
		final var fluidBlockRhs = Fluid.assembleBlockRhs(fluidSystem, getDT());
		blocks.put(new IntCoordinates(0, 0), fluidBlockRhs._1);
		rhs.addSmallVectorAt(fluidBlockRhs._2, 0);
		offset += fluidBlockRhs._1.getCols();
		for (int i = 0; i < getParticles().size(); i++)
		{
			offset = addParticleBlocks(blocks, rhs, particleSystems, offset, i);
		}
		final BlockSparseMatrix systemMatrix = new BlockSparseMatrix(blocks, rhs.getLength(), rhs.getLength());
		final Tuple2<BlockSparseMatrix, DenseVector> linearSystem = applyBoundaryValues(systemMatrix,
		                                                                                rhs,
		                                                                                time);
		return new Tuple4<>(fluidSystem, particleSystems, linearSystem._1, linearSystem._2);
	}
	
	default int getTotalSystemSize()
	{
		return getFluid().getSystemSize() + getParticles().stream()
		                                                  .mapToInt(p -> p.getSystemSize() + p.getLagrangeSize())
		                                                  .sum();
	}
	
	@NotNull
	default Tuple2<FluidIterate, List<ParticleIterate>> createNewIterates(final List<ParticleIterate> previousParticleStates,
	                                                                      final Vector solution)
	{
		int offset;
		final FluidIterate ret = new FluidIterate(solution.slice(0, getFluid().getSystemSize()));
		offset = getFluid().getSystemSize();
		final List<ParticleIterate> iterates = new ArrayList<>();
		for (int i = 0; i < getParticles().size(); i++)
		{
			final ParticleIterate it
				= new ParticleIterate(previousParticleStates.get(i),
				                      solution.slice(offset,
				                                     offset + getParticles().get(i)
				                                                            .getSystemSize())
				                              .mul(getDT()),
				                      solution.slice(offset + getParticles().get(i)
				                                                            .getSystemSize(),
				                                     offset + getParticles().get(i)
				                                                            .getSystemSize()
					                                     + getParticles().get(i)
					                                                     .getLagrangeSize()));
			iterates.add(it);
			offset += getParticles().get(i)
			                        .getSystemSize() + getParticles().get(i)
			                                                         .getLagrangeSize();
		}
		return new Tuple2<>(ret, iterates);
	}
	
	private Tuple2<BlockSparseMatrix, DenseVector> applyBoundaryValues(final BlockSparseMatrix systemMatrix,
	                                                                   final DenseVector rhs, final double t)
	{
		final SparseMatrix ret = systemMatrix.toSparse();
		final DenseVector retRhs = new DenseVector(rhs);
		final Int2DoubleMap nodeValues = new Int2DoubleArrayMap();
		final Int2DoubleMap fluidDirichletNodeValues = getFluid().getDirichletNodeValues(t);
		nodeValues.putAll(fluidDirichletNodeValues);
		for (int i = 0; i < getParticles().size(); i++)
		{
			final Int2DoubleMap particleDirichletNodeValues = getParticles().get(i)
			                                                                .getDirichletNodeValues(t);
			final int finalI = i;
			final int particleSpaceSize = getParticles().get(i)
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
		final var particleBlockRhs = getParticles().get(i)
		                                           .getBlockRhs(particleSystems.get(i), getDT());
		blocks.put(new IntCoordinates(offset, offset), particleBlockRhs._1);
		rhs.addSmallVectorAt(particleBlockRhs._2, offset);
		final var particleBackgroundLagrangeBlock =
			getParticles().get(i)
			              .getLagrangeBackgroundBlock(particleSystems.get(i), getFluid(), getDT());
		blocks.put(new IntCoordinates(0, offset), particleBackgroundLagrangeBlock);
		final var particleBackgroundLagrangeBlockTranspose =
			getParticles().get(i)
			              .getLagrangeBackgroundBlockTranspose(particleSystems.get(i),
			                                                   getFluid(),
			                                                   getDT());
		blocks.put(new IntCoordinates(offset, 0), particleBackgroundLagrangeBlockTranspose);
		offset += particleBlockRhs._1.getCols();
		return offset;
	}
	
	void postIterationCallback(final FluidIterate fluidState,
	                           final List<ParticleIterate> particleStates,
	                           final double time);
	
	void show(final FluidIterate fluidState,
	          final List<ParticleIterate> particleStates,
	          int iteration);
}
