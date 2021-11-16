package dlm;

import basic.*;
import distorted.*;
import linalg.*;
import mixed.*;
import tensorproduct.CartesianGridSpace;
import tensorproduct.TPVectorCellIntegral;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.function.Function;
import java.util.function.ToDoubleFunction;

public abstract class HyperbolicCartesianDistorted<ST extends MixedShapeFunction<TPCell, TPFace, PF, VF>,
	PF extends ScalarShapeFunction<TPCell, TPFace>, VF extends VectorShapeFunction<TPCell, TPFace>>
	extends FullyImplicitHyperbolicIntegrator
{
	final int dimension;
	final CartesianGridSpace<ST, MixedValue, MixedGradient, MixedHessian> backgroundSpace;
	final List<DistortedVectorSpace> particleSpaces;
	final BlockSparseMatrix twoDerivativeMatrix;
	final BlockSparseMatrix oneDerivativeMatrix;
	final BlockSparseMatrix zeroDerivativeMatrix;
	final List<Double> rhoS;
	final double rhoF;
	
	public HyperbolicCartesianDistorted(final double dt,
	                                    final int timeSteps,
	                                    final CartesianGridSpace<ST, MixedValue, MixedGradient, MixedHessian> backgroundSpace,
	                                    final List<DistortedVectorSpace> particleSpaces,
	                                    final List<Double> particleDensity, final double backgroundDensity)
	{
		super(dt, timeSteps);
		this.backgroundSpace = backgroundSpace;
		this.particleSpaces = particleSpaces;
		this.dimension = backgroundSpace.getDimension();
		this.rhoS = particleDensity;
		this.rhoF = backgroundDensity;
		final Map<IntCoordinates, SparseMatrix> twoDerivativeBlocks = new HashMap<>();
		final Map<IntCoordinates, SparseMatrix> oneDerivativeBlocks = new HashMap<>();
		final Map<IntCoordinates, SparseMatrix> zeroDerivativeBlocks = new HashMap<>();
		putBackgroundBlocks(twoDerivativeBlocks, oneDerivativeBlocks, zeroDerivativeBlocks);
		putParticleBlocks(particleSpaces, twoDerivativeBlocks, oneDerivativeBlocks, zeroDerivativeBlocks);
		twoDerivativeMatrix = new BlockSparseMatrix(twoDerivativeBlocks);
		oneDerivativeMatrix = new BlockSparseMatrix(oneDerivativeBlocks);
		zeroDerivativeMatrix = new BlockSparseMatrix(zeroDerivativeBlocks);
	}
	
	private void putParticleBlocks(final List<DistortedVectorSpace> particleSpaces,
	                               final Map<IntCoordinates, SparseMatrix> twoDerivativeBlocks,
	                               final Map<IntCoordinates, SparseMatrix> oneDerivativeBlocks,
	                               final Map<IntCoordinates, SparseMatrix> zeroDerivativeBlocks)
	{
		int offset = backgroundSpace.getShapeFunctions()
		                            .size();
		for (int i = 0; i < particleSpaces.size(); i++)
		{
			putParticleBlocks(twoDerivativeBlocks, oneDerivativeBlocks, zeroDerivativeBlocks, offset, i);
			offset += particleSpaces.get(i)
			                        .getShapeFunctions()
			                        .size() * 2;
		}
	}
	
	private void putParticleBlocks(final Map<IntCoordinates, SparseMatrix> twoDerivativeBlocks,
	                               final Map<IntCoordinates, SparseMatrix> oneDerivativeBlocks,
	                               final Map<IntCoordinates, SparseMatrix> zeroDerivativeBlocks,
	                               final int offset,
	                               final int particleId)
	{
		final DistortedVectorSpace particleSpace = particleSpaces.get(particleId);
		final int nParticle = particleSpace.getShapeFunctions()
		                                   .size();
		final int nBackgroundVelocities = getnBackgroundVelocities();
		final SparseMatrix Cs = getCs(particleId);
		twoDerivativeBlocks.put(new IntCoordinates(offset, offset), getMs(particleId));
		twoDerivativeBlocks.put(new IntCoordinates(offset + nParticle, offset + nParticle),
		                        new SparseMatrix(nParticle,
		                                         nParticle));
		oneDerivativeBlocks.put(new IntCoordinates(offset + nParticle, offset), Cs);
		oneDerivativeBlocks.put(new IntCoordinates(offset, offset), new SparseMatrix(nParticle, nParticle));
		oneDerivativeBlocks.put(new IntCoordinates(offset + nParticle, offset + nParticle),
		                        new SparseMatrix(nParticle,
		                                         nParticle));
		zeroDerivativeBlocks.put(new IntCoordinates(offset, offset + nParticle), Cs.transpose());
		
		zeroDerivativeBlocks.put(new IntCoordinates(offset, offset), getAs());
		zeroDerivativeBlocks.put(new IntCoordinates(offset + nParticle, offset + nParticle),
		                         new SparseMatrix(nParticle,
		                                          nParticle));
		zeroDerivativeBlocks.put(new IntCoordinates(0, offset + nParticle),
		                         new SparseMatrix(nBackgroundVelocities, nParticle));
		zeroDerivativeBlocks.put(new IntCoordinates(offset + nParticle, 0),
		                         new SparseMatrix(nParticle, nBackgroundVelocities));
	}
	
	private void putBackgroundBlocks(final Map<IntCoordinates, SparseMatrix> twoDerivativeBlocks,
	                                 final Map<IntCoordinates, SparseMatrix> oneDerivativeBlocks,
	                                 final Map<IntCoordinates, SparseMatrix> zeroDerivativeBlocks)
	{
		final int nBackgroundVelocities = getnBackgroundVelocities();
		final int nBackgroundPressures = getnBackgroundPressures();
		twoDerivativeBlocks.put(new IntCoordinates(0, 0),
		                        new SparseMatrix(nBackgroundVelocities, nBackgroundVelocities));
		twoDerivativeBlocks.put(new IntCoordinates(nBackgroundVelocities, nBackgroundVelocities),
		                        new SparseMatrix(nBackgroundPressures, nBackgroundPressures));
		oneDerivativeBlocks.put(new IntCoordinates(0, 0), getMf());
		oneDerivativeBlocks.put(new IntCoordinates(nBackgroundVelocities, nBackgroundVelocities),
		                        new SparseMatrix(nBackgroundPressures, nBackgroundPressures));
		final SparseMatrix ABf = getABf();
		final SparseMatrix Af = ABf.slice(new IntCoordinates(0, 0), new IntCoordinates(nBackgroundVelocities,
		                                                                               nBackgroundVelocities));
		final SparseMatrix Bf = ABf.slice(new IntCoordinates(nBackgroundVelocities, 0),
		                                  new IntCoordinates(nBackgroundPressures,
		                                                     nBackgroundVelocities));
		zeroDerivativeBlocks.put(new IntCoordinates(0, 0),
		                         Af);
		zeroDerivativeBlocks.put(new IntCoordinates(nBackgroundVelocities, 0),
		                         Bf);
		zeroDerivativeBlocks.put(new IntCoordinates(0, nBackgroundVelocities),
		                         Bf.transpose());
		zeroDerivativeBlocks.put(new IntCoordinates(nBackgroundVelocities, nBackgroundVelocities),
		                         new SparseMatrix(nBackgroundPressures, nBackgroundPressures));
	}
	
	protected abstract SparseMatrix getABf();
	
	protected abstract SparseMatrix getAs();
	
	private SparseMatrix getCs(final int particleId)
	{
		final DistortedVectorCellIntegral particleH1
			= new DistortedVectorCellIntegral(DistortedVectorCellIntegral.H1);
		final SparseMatrix s = new SparseMatrix(particleSpaces.get(particleId)
		                                                      .getShapeFunctions()
		                                                      .size(),
		                                        particleSpaces.get(particleId)
		                                                      .getShapeFunctions()
		                                                      .size());
		particleSpaces.get(particleId)
		              .writeCellIntegralsToMatrix(List.of(particleH1), s);
		return s;
	}
	
	private SparseMatrix getMf()
	{
		final TPVectorCellIntegral<VF> fluidMass = new TPVectorCellIntegral<>(rhoF,
		                                                                      TPVectorCellIntegral.VALUE_VALUE);
		final MixedCellIntegral<TPCell, PF, VF, ST> mixedFluidMass
			= MixedCellIntegral.fromVelocityIntegral(fluidMass);
		final SparseMatrix s = new SparseMatrix(backgroundSpace.getShapeFunctions()
		                                                       .size(),
		                                        backgroundSpace.getShapeFunctions()
		                                                       .size());
		backgroundSpace.writeCellIntegralsToMatrix(List.of(mixedFluidMass), s);
		return s.slice(new IntCoordinates(0, 0), new IntCoordinates(getnBackgroundVelocities(),
		                                                            getnBackgroundVelocities()));
	}
	
	private SparseMatrix getMs(final int particleId)
	{
		final DistortedVectorCellIntegral particleMass = new DistortedVectorCellIntegral(rhoS.get(particleId),
		                                                                                 DistortedVectorCellIntegral.VALUE_VALUE);
		final SparseMatrix s = new SparseMatrix(particleSpaces.get(particleId)
		                                                      .getShapeFunctions()
		                                                      .size(),
		                                        particleSpaces.get(particleId)
		                                                      .getShapeFunctions()
		                                                      .size());
		particleSpaces.get(particleId)
		              .writeCellIntegralsToMatrix(List.of(particleMass), s);
		return s;
	}
	
	private SparseMatrix getCf(final int particleId)
	{
		final SparseMatrix Cf = new SparseMatrix(zeroDerivativeMatrix.getBlockSizes()[particleId * 2],
		                                         zeroDerivativeMatrix.getBlockSizes()[0]);
		backgroundSpace.getShapeFunctions()
		               .entrySet()
		               .stream()
		               .parallel()
		               .filter(entry -> entry.getValue()
		                                     .hasVelocityFunction())
		               .forEach(entry ->
		                        {
			                        final int index = entry.getKey();
			                        final VF velocityShapeFunction =
				                        entry.getValue()
				                             .getVelocityShapeFunction();
			                        final DistortedVectorFunction concatenation =
				                        DistortedVectorFunctionOnCells.concatenate(velocityShapeFunction,
				                                                                   getDisplacement(
					                                                                   particleId));
			                        final DistortedVectorDistortedRightHandSideIntegral backGroundH1 =
				                        new DistortedVectorDistortedRightHandSideIntegral
					                        (concatenation,
					                         DistortedVectorDistortedRightHandSideIntegral.H1);
			                        final DenseVector column =
				                        new DenseVector(Cf.getRows());
			                        particleSpaces.get(particleId)
			                                      .writeCellIntegralsToRhs(List.of(backGroundH1), column);
			                        Cf.addColumn(column, index);
		                        });
		return Cf;
	}
	
	private DistortedVectorFunctionOnCells getDisplacement(final int particleId)
	{
		
		final DenseVector coefficients
			= currentIterate.slice(zeroDerivativeMatrix.getBlockStarts()[particleId * 2],
			                       zeroDerivativeMatrix.getBlockEnds()[particleId * 2]);
		return new DistortedVectorFESpaceFunction(particleSpaces.get(particleId)
		                                                        .getShapeFunctions(),
		                                          coefficients);
	}
	
	private int getnBackgroundPressures()
	{
		return (int) backgroundSpace.getShapeFunctions()
		                            .values()
		                            .stream()
		                            .filter(MixedFunction::hasPressureFunction)
		                            .count();
	}
	
	private int getnBackgroundVelocities()
	{
		return (int) backgroundSpace.getShapeFunctions()
		                            .values()
		                            .stream()
		                            .filter(MixedFunction::hasVelocityFunction)
		                            .count();
	}
	
	abstract protected Function<CoordinateVector, CoordinateVector> getInitialFluidVelocity();
	
	abstract protected ToDoubleFunction<CoordinateVector> getInitialPressure();
	
	abstract protected Function<CoordinateVector, CoordinateVector> getInitialParticleVelocity(int particleId);
	
	@Override
	protected DenseVector initializeInitialIterate()
	{
		final DenseVector ret = new DenseVector(zeroDerivativeMatrix.getRows());
		addInitialUp(ret);
		for (int i = 0; i < particleSpaces.size(); i++)
		{
			addInitialParticleDisplacement(ret, i);
		}
		return ret;
	}
	
	private void addInitialParticleDisplacement(final DenseVector ret, final int particleId)
	{
		final VectorFunction identity = VectorFunction.fromLambda(x -> x, dimension, dimension);
		final DenseVector initialDisplacement
			= new DenseVector(zeroDerivativeMatrix.getBlockSizes()[particleId + 2]);
		particleSpaces.get(particleId)
		              .getShapeFunctions()
		              .forEach((key, func) -> initialDisplacement.add(func.getNodeFunctional()
		                                                                  .evaluate(identity), key));
		ret.addSmallVectorAt(initialDisplacement, zeroDerivativeMatrix.getBlockStarts()[particleId + 2]);
	}
	
	private void addInitialUp(final DenseVector ret)
	{
		final MixedFunction initialUpFunction
			= new MixedFunction(ScalarFunction.fromLambda(getInitialPressure(),
			                                              dimension),
			                    VectorFunction.fromLambda(getInitialFluidVelocity(),
			                                              dimension,
			                                              dimension));
		final DenseVector initialUp = new DenseVector(getnBackgroundVelocities());
		backgroundSpace.getShapeFunctions()
		               .forEach((key, value) ->
			                        initialUp.add(value.getNodeFunctional()
			                                           .evaluate(initialUpFunction), key));
		ret.addSmallVectorAt(initialUp, 0);
	}
	
	@Override
	protected DenseVector initializeInitialDerivative()
	{
		final DenseVector ret = new DenseVector(zeroDerivativeMatrix.getRows());
		for (int i = 0; i < particleSpaces.size(); i++)
		{
			addInitialParticleVelocity(ret, i);
		}
		return ret;
	}
	
	private void addInitialParticleVelocity(final DenseVector ret, final int particleId)
	{
		final DenseVector initialVelocity
			= new DenseVector(zeroDerivativeMatrix.getBlockSizes()[particleId + 2]);
		particleSpaces.get(particleId)
		              .getShapeFunctions()
		              .forEach((key, func) ->
			                       initialVelocity.add(func.getNodeFunctional()
			                                               .evaluate(VectorFunction.fromLambda(
				                                               getInitialParticleVelocity(particleId),
				                                               dimension,
				                                               dimension)), key));
		ret.addSmallVectorAt(initialVelocity, zeroDerivativeMatrix.getBlockStarts()[particleId + 2]);
	}
	
	@Override
	protected VectorMultiplyable getDoubleDerivativeOperator()
	{
		return twoDerivativeMatrix;
	}
	
	@Override
	protected VectorMultiplyable getSingleDerivativeOperator()
	{
		return oneDerivativeMatrix;
	}
	
	@Override
	protected VectorMultiplyable getNoDerivativeOperator()
	{
		final SparseMatrix As = zeroDerivativeMatrix.getBlockMatrix(0, 0);
		
		for (int i = 0; i < particleSpaces.size(); i++)
		{
			final SparseMatrix Cf = getCf(i);
			zeroDerivativeMatrix.getBlockMatrix(i, 0)
			                    .overrideBy(Cf);
			zeroDerivativeMatrix.getBlockMatrix(0, i)
			                    .overrideBy(Cf.transpose());
		}
		return zeroDerivativeMatrix;
	}
}
