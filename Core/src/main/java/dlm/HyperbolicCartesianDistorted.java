package dlm;

import basic.*;
import com.google.common.collect.ImmutableList;
import distorted.DistortedGridSpace;
import distorted.DistortedVectorFESpaceFunction;
import distorted.DistortedVectorFunctionOnCells;
import distorted.DistortedVectorShapeFunction;
import distorted.geometry.DistortedCell;
import linalg.Vector;
import linalg.*;
import mixed.*;
import scala.Function2;
import tensorproduct.CartesianGridSpace;
import tensorproduct.ContinuousTPVectorFunction;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

import java.util.*;
import java.util.function.Function;
import java.util.function.Predicate;
import java.util.function.ToDoubleFunction;

public abstract class HyperbolicCartesianDistorted<ST extends ComposeMixedShapeFunction<TPCell, TPFace, PF, VF>,
	PF extends ScalarShapeFunction<TPCell, TPFace>, VF extends VectorShapeFunction<TPCell, TPFace>>
	extends FullyImplicitHyperbolicIntegrator
{
	final public int dimension;
	final public CartesianGridSpace<ST, MixedValue, MixedGradient, MixedHessian> backgroundSpace;
	final public ImmutableList<DistortedGridSpace<DistortedVectorShapeFunction, CoordinateVector, CoordinateMatrix, CoordinateTensor>>
		particleSpaces;
	final public BlockSparseMatrix twoDerivativeMatrix;
	final public BlockSparseMatrix oneDerivativeMatrix;
	final public BlockSparseMatrix zeroDerivativeMatrix;
	final int[] blockStarts;
	SparseMatrix fixedBackGroundVelocityMatrix;
	final MetricWindow metricWindow = MetricWindow.getInstance();
	final CountMetric cm;
	
	public HyperbolicCartesianDistorted(final double dt,
	                                    final int timeSteps,
	                                    final CartesianGridSpace<ST, MixedValue, MixedGradient, MixedHessian> backgroundSpace,
	                                    final List<DistortedGridSpace<DistortedVectorShapeFunction, CoordinateVector, CoordinateMatrix, CoordinateTensor>> particleSpaces)
	{
		super(dt, timeSteps);
		this.backgroundSpace = backgroundSpace;
		this.particleSpaces = ImmutableList.copyOf(particleSpaces);
		this.dimension = backgroundSpace.getDimension();
		final Map<IntCoordinates, SparseMatrix> twoDerivativeBlocks = new HashMap<>();
		final Map<IntCoordinates, SparseMatrix> oneDerivativeBlocks = new HashMap<>();
		final Map<IntCoordinates, SparseMatrix> zeroDerivativeBlocks = new HashMap<>();
		blockStarts = new int[2 * particleSpaces.size() + 3];
		blockStarts[0] = 0;
		blockStarts[1] = getnBackgroundVelocities();
		blockStarts[2] = getnBackgroundVelocities() + getnBackgroundPressures();
		for (int i = 3; i < 2 + 2 * particleSpaces.size(); i++)
		{
			final int particleId = (i - 3) / 2;
			blockStarts[i] = blockStarts[i - 1] + particleSpaces.get(particleId)
			                                                    .getShapeFunctions()
			                                                    .size();
		}
		blockStarts[2 * particleSpaces.size() + 2] = getSystemSize();
		System.out.println(Arrays.toString(blockStarts));
		addBackgroundBlocks(twoDerivativeBlocks, oneDerivativeBlocks, zeroDerivativeBlocks);
		addAllParticleBlocks(twoDerivativeBlocks, oneDerivativeBlocks, zeroDerivativeBlocks);
		final int size = getSystemSize();
		twoDerivativeMatrix = new BlockSparseMatrix(twoDerivativeBlocks, size, size);
		oneDerivativeMatrix = new BlockSparseMatrix(oneDerivativeBlocks, size, size);
		zeroDerivativeMatrix = new BlockSparseMatrix(zeroDerivativeBlocks, size, size);
		cm = new CountMetric(getnBackgroundVelocities());
		metricWindow.addMetric(cm);
	}
	
	protected int getSystemSize()
	{
		return backgroundSpace.getShapeFunctions()
		                      .size() + particleSpaces.stream()
		                                              .mapToInt(sp -> sp.getShapeFunctions()
		                                                                .size())
		                                              .sum() * 2;
	}
	
	private void addAllParticleBlocks(final Map<IntCoordinates, SparseMatrix> twoDerivativeBlocks,
	                                  final Map<IntCoordinates, SparseMatrix> oneDerivativeBlocks,
	                                  final Map<IntCoordinates, SparseMatrix> zeroDerivativeBlocks)
	{
		for (int particleId = 0; particleId < particleSpaces.size(); particleId++)
		{
			addSingleParticleBlocks(twoDerivativeBlocks,
			                        oneDerivativeBlocks,
			                        zeroDerivativeBlocks,
			                        particleId);
		}
	}
	
	private void addBackgroundBlocks(final Map<IntCoordinates, SparseMatrix> twoDerivativeBlocks,
	                                 final Map<IntCoordinates, SparseMatrix> oneDerivativeBlocks,
	                                 final Map<IntCoordinates, SparseMatrix> zeroDerivativeBlocks)
	{
		final int nBackgroundVelocities = getnBackgroundVelocities();
		final int nBackgroundPressures = getnBackgroundPressures();
		final SparseMatrix mixedBackgroundMassMatrix = evaluateBackgroundMassIntegrals();
		final SparseMatrix velocityMassMatrix = mixedBackgroundMassMatrix.slice(new IntCoordinates(0, 0),
		                                                                        new IntCoordinates(
			                                                                        nBackgroundVelocities,
			                                                                        nBackgroundVelocities));
		oneDerivativeBlocks.put(new IntCoordinates(0, 0), velocityMassMatrix);
		final SparseMatrix mixedBackgroundMatrix = evaluateBackgroundIntegrals();
		fixedBackGroundVelocityMatrix = mixedBackgroundMatrix.slice(new IntCoordinates(0, 0),
		                                                            new IntCoordinates(nBackgroundVelocities,
		                                                                               nBackgroundVelocities));
		final SparseMatrix incompressibilityMatrix = mixedBackgroundMatrix.slice(new IntCoordinates(
			                                                                         nBackgroundVelocities,
			                                                                         0),
		                                                                         new IntCoordinates(
			                                                                         nBackgroundVelocities + nBackgroundPressures,
			                                                                         nBackgroundVelocities));
		zeroDerivativeBlocks.put(new IntCoordinates(0, 0), fixedBackGroundVelocityMatrix);
		zeroDerivativeBlocks.put(new IntCoordinates(nBackgroundVelocities, 0),
		                         incompressibilityMatrix);
		zeroDerivativeBlocks.put(new IntCoordinates(0, nBackgroundVelocities),
		                         incompressibilityMatrix.transpose());
	}
	
	private void addSingleParticleBlocks(final Map<IntCoordinates, SparseMatrix> twoDerivativeBlocks,
	                                     final Map<IntCoordinates, SparseMatrix> oneDerivativeBlocks,
	                                     final Map<IntCoordinates, SparseMatrix> zeroDerivativeBlocks,
	                                     final int particleId)
	{
		final DistortedGridSpace<DistortedVectorShapeFunction, CoordinateVector, CoordinateMatrix, CoordinateTensor>
			particleSpace = particleSpaces.get(particleId);
		final int particleBlockSize = particleSpace.getShapeFunctions()
		                                           .size();
		final int nBackgroundVelocities = getnBackgroundVelocities();
		final SparseMatrix particleMass = evaluateParticleMassIntegrals(particleId);
		final SparseMatrix elasticityMatrix = evaluateParticleIntegrals(particleId);
		final SparseMatrix particleLagrange = evaluateParticleLagrangeIntegrals(particleId);
		final SparseMatrix backgroundLagrangePlaceHolder = new SparseMatrix(nBackgroundVelocities,
		                                                                    particleBlockSize);
		final int block = 2 * particleId + 2;
		twoDerivativeBlocks.put(new IntCoordinates(blockStarts[block], blockStarts[block]),
		                        particleMass);
		oneDerivativeBlocks.put(new IntCoordinates(blockStarts[block + 1], blockStarts[block]),
		                        particleLagrange.transpose());
		zeroDerivativeBlocks.put(new IntCoordinates(blockStarts[block], blockStarts[block + 1]),
		                         particleLagrange);
		zeroDerivativeBlocks.put(new IntCoordinates(blockStarts[block], blockStarts[block]),
		                         elasticityMatrix);
		zeroDerivativeBlocks.put(new IntCoordinates(0, blockStarts[block + 1]),
		                         backgroundLagrangePlaceHolder);
		zeroDerivativeBlocks.put(new IntCoordinates(blockStarts[block + 1], 0),
		                         backgroundLagrangePlaceHolder.transpose());
	}
	
	protected abstract List<CellIntegral<TPCell, ST>> getBackgroundIntegrals();
	
	private SparseMatrix evaluateBackgroundIntegrals()
	{
		final int n = blockStarts[2] - blockStarts[0];
		final SparseMatrix ret = new SparseMatrix(n, n);
		backgroundSpace.writeCellIntegralsToMatrix(getBackgroundIntegrals(), ret);
		return ret;
	}
	
	protected List<CellIntegral<TPCell, ST>> getSemiImplicitBackgroundIntegrals(final VectorFunctionOnCells<TPCell,
		TPFace> u)
	{
		return new ArrayList<>();
	}
	
	protected SparseMatrix evaluateSemiImplicitBackGroundIntegrals()
	{
		final int n = blockStarts[2] - blockStarts[0];
		final SparseMatrix ret = new SparseMatrix(n, n);
		backgroundSpace.writeCellIntegralsToMatrix(getSemiImplicitBackgroundIntegrals(getUp().getVelocityFunction()),
		                                           ret);
		return ret;
	}
	
	protected abstract List<CellIntegral<DistortedCell, DistortedVectorShapeFunction>> getParticleIntegrals(int particleId);
	
	private SparseMatrix evaluateParticleIntegrals(final int particleId)
	{
		final int n = blockStarts[2 * particleId + 2 + 1] - blockStarts[2 * particleId + 2];
		final SparseMatrix ret = new SparseMatrix(n, n);
		particleSpaces.get(particleId)
		              .writeCellIntegralsToMatrix(getParticleIntegrals(particleId), ret);
		return ret;
	}
	
	protected abstract List<CellIntegral<DistortedCell, DistortedVectorShapeFunction>> getParticleLagrangeIntegrals(
		int particleId);
	
	private SparseMatrix evaluateParticleLagrangeIntegrals(final int particleId)
	{
		final int n = blockStarts[2 * particleId + 2 + 1] - blockStarts[2 * particleId + 2];
		final SparseMatrix ret = new SparseMatrix(n, n);
		particleSpaces.get(particleId)
		              .writeCellIntegralsToMatrix(getParticleLagrangeIntegrals(particleId), ret);
		return ret;
	}
	
	protected abstract List<RightHandSideIntegral<DistortedCell, DistortedVectorShapeFunction>> getBackgroundLagrangeIntegrals(
		ST function,
		DistortedVectorFunctionOnCells X,
		int particleId);
	
	private SparseMatrix evaluateBackgroundLagrangeIntegrals(final int particleId)
	{
		final int nCols = blockStarts[2 * particleId + 4] - blockStarts[2 * particleId + 3];
		final int nRows = blockStarts[1] - blockStarts[0];
		final SparseMatrix ret = new SparseMatrix(nRows, nCols);
		final DistortedVectorFunctionOnCells X = getX(particleId);
		final int nEulerCells = (int) (Math.sqrt(backgroundSpace.getCells()
		                                                        .size()));
		final double diam = particleSpaces.get(particleId)
		                                  .getMaxDiam();
		
		backgroundSpace.getShapeFunctions()
		               .entrySet()
		               .stream()
		               .parallel()
		               .filter(e -> e.getValue()
		                             .hasVelocityFunction())
		               .forEach(entry ->
		                        {
			                        cm.increment();
			                        final ST backGroundFunction = entry.getValue();
			                        final int rowIndex = entry.getKey();
			                        final DenseVector row = new DenseVector(nCols);
			                        particleSpaces
				                        .get(particleId)
				                        .writeCellIntegralsToRhs(
					                        getBackgroundLagrangeIntegrals(
						                        backGroundFunction,
						                        X,
						                        particleId),
					                        row,
					                        (cell, fun)
						                        -> ((ContinuousTPVectorFunction) backGroundFunction
						                        .getVelocityFunction()).getNodeFunctionalPoint()
						                                               .sub(cell.center())
						                                               .euclidianNorm() < 2. / nEulerCells + diam);
			                        ret.addRow(row, rowIndex);
		                        });
		cm.reset();
		return ret;
	}
	
	public DistortedVectorFunctionOnCells getX(final int particleId)
	{
		final DistortedVectorFESpaceFunction displacement
			= new DistortedVectorFESpaceFunction(particleSpaces.get(
				                                                   particleId)
			                                                   .getShapeFunctions(),
			                                     getParticleIterate(
				                                     currentIterate,
				                                     particleId));
		return new DistortedVectorFunctionOnCells()
		{
			@Override
			public int getRangeDimension()
			{
				return 2;
			}
			
			@Override
			public CoordinateVector valueOnReferenceCell(final CoordinateVector pos,
			                                             final DistortedCell cell)
			{
				return cell.transformFromReferenceCell(pos)
				           .add(displacement.valueOnReferenceCell(pos
					           , cell));
			}
			
			@Override
			public CoordinateMatrix gradientOnReferenceCell(final CoordinateVector pos,
			                                                final DistortedCell cell)
			{
				return CoordinateDenseMatrix.identity(2)
				                            .add(displacement.gradientOnReferenceCell(pos
					                            , cell));
			}
			
			@Override
			public int getDomainDimension()
			{
				return 2;
			}
			
			@Override
			public CoordinateVector value(final CoordinateVector pos)
			{
				throw new UnsupportedOperationException("not implemented yet");
			}
		};
	}
	
	public MixedFunctionOnCells<TPCell, TPFace> getUp()
	{
		return new MixedTPFESpaceFunction<>(backgroundSpace.getShapeFunctions(),
		                                    getBackGroundIterate(currentIterate));
	}
	
	protected abstract List<CellIntegral<TPCell, ST>> getBackgroundMassIntegrals();
	
	private SparseMatrix evaluateBackgroundMassIntegrals()
	{
		final int n = blockStarts[2] - blockStarts[0];
		final SparseMatrix ret = new SparseMatrix(n, n);
		backgroundSpace.writeCellIntegralsToMatrix(getBackgroundMassIntegrals(), ret);
		return ret;
	}
	
	protected abstract List<CellIntegral<DistortedCell, DistortedVectorShapeFunction>> getParticleMassIntegrals(int particleId);
	
	private SparseMatrix evaluateParticleMassIntegrals(final int particleId)
	{
		final int n = blockStarts[2 * particleId + 2 + 1] - blockStarts[2 * particleId + 2];
		final SparseMatrix ret = new SparseMatrix(n, n);
		particleSpaces.get(particleId)
		              .writeCellIntegralsToMatrix(getParticleMassIntegrals(particleId), ret);
		return ret;
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
		final DenseVector ret = new DenseVector(getSystemSize());
		addInitialUp(ret);
		for (int i = 0; i < particleSpaces.size(); i++)
		{
			addInitialParticleDisplacement(ret, i);
		}
		return ret;
	}
	
	private void addInitialParticleDisplacement(final DenseVector ret, final int particleId)
	{
		final VectorFunction identity = VectorFunction.fromLambda(x -> x.mul(0), dimension, dimension);
		final DenseVector initialDisplacement
			= new DenseVector(zeroDerivativeMatrix.getBlockSizes()[2 * particleId + 2]);
		particleSpaces.get(particleId)
		              .getShapeFunctions()
		              .forEach((key, func) -> initialDisplacement.add(func.getNodeFunctional()
		                                                                  .evaluate(identity), key));
		ret.addSmallVectorAt(initialDisplacement, zeroDerivativeMatrix.getBlockStarts()[particleId + 2]);
	}
	
	private void addInitialUp(final DenseVector ret)
	{
		final MixedFunction initialUpFunction
			= new ComposedMixedFunction(ScalarFunction.fromLambda(getInitialPressure(),
			                                                      dimension),
			                            VectorFunction.fromLambda(getInitialFluidVelocity(),
			                                                      dimension,
			                                                      dimension));
		final DenseVector initialUp = new DenseVector(getnBackgroundVelocities() + getnBackgroundPressures());
		backgroundSpace.getShapeFunctions()
		               .forEach((key, value) ->
			                        initialUp.add(value.getNodeFunctional()
			                                           .evaluate(initialUpFunction), key));
		ret.addSmallVectorAt(initialUp, 0);
	}
	
	public DenseVector getBackGroundIterate(final Vector vector)
	{
		return vector.slice(0, blockStarts[2]);
	}
	
	public DenseVector getBackGroundVelocityIterate(final Vector vector)
	{
		return vector.slice(0, blockStarts[1]);
	}
	
	public DenseVector getBackGroundPressureIterate(final Vector vector)
	{
		return vector.slice(blockStarts[1], blockStarts[2]);
	}
	
	public DenseVector getParticleIterate(final Vector vector, final int particleId)
	{
		return vector.slice(blockStarts[2 * particleId + 2], blockStarts[2 * particleId + 3]);
	}
	
	public DenseVector getAllParticleIterate(final Vector vector, final int particleId)
	{
		return vector.slice(blockStarts[2 * particleId + 2], blockStarts[2 * particleId + 4]);
	}
	
	public DenseVector getParticleTransferIterate(final Vector vector, final int particleId)
	{
		return vector.slice(blockStarts[2 * particleId + 3], blockStarts[2 * particleId + 4]);
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
			= new DenseVector(zeroDerivativeMatrix.getBlockSizes()[2 * particleId + 2]);
		final VectorFunction initialParticleVelocity = VectorFunction.fromLambda(
			getInitialParticleVelocity(particleId), dimension, dimension);
		particleSpaces.get(particleId)
		              .getShapeFunctions()
		              .forEach((key, func) ->
			                       initialVelocity.add(func.getNodeFunctional()
			                                               .evaluate(initialParticleVelocity), key));
		ret.addSmallVectorAt(initialVelocity, zeroDerivativeMatrix.getBlockStarts()[2 * particleId + 2]);
	}
	
	@Override
	protected Matrix getDoubleDerivativeOperator()
	{
		return twoDerivativeMatrix;
	}
	
	@Override
	protected Matrix getSingleDerivativeOperator()
	{
		return oneDerivativeMatrix;
	}
	
	@Override
	protected Matrix getNoDerivativeOperator()
	{
		final SparseMatrix As =
			evaluateSemiImplicitBackGroundIntegrals().slice(new IntCoordinates(0, 0),
			                                                new IntCoordinates(getnBackgroundVelocities(),
			                                                                   getnBackgroundVelocities()))
			                                         .add(fixedBackGroundVelocityMatrix);
		zeroDerivativeMatrix.getBlockMatrix(0, 0)
		                    .overrideBy(As);
		
		for (int i = 0; i < particleSpaces.size(); i++)
		{
			final SparseMatrix Cf = evaluateBackgroundLagrangeIntegrals(i);
			zeroDerivativeMatrix.getBlocks()
			                    .get(new IntCoordinates(blockStarts[2 * i + 3], blockStarts[0]))
			                    .overrideBy(Cf.transpose());
			zeroDerivativeMatrix.getBlocks()
			                    .get(new IntCoordinates(blockStarts[0], blockStarts[2 * i + 3]))
			                    .overrideBy(Cf);
		}
		return zeroDerivativeMatrix;
	}
	
	protected Matrix getSchurComplement()
	{
		final SparseMatrix ret = new SparseMatrix(getAllBackgroundMatrices());
		for (int i = 0; i < particleSpaces.size(); i++)
		{
			final DenseMatrix D = new DenseMatrix(getAllParticleMatrices(i));
			final Matrix B = getAllBackgroundTransferMatrices(i);
			final Matrix C = getAllBackgroundTransferTransposeMatrices(i);
			final Matrix Dinv = D.inverse();
			final Matrix BD = B.mmMul(Dinv);
			final Matrix localMatrix = BD.mmMul(C);
			ret.subInPlace(localMatrix);
		}
		return ret;
	}
	
	protected Matrix getAllParticleMatrices(final int particleId)
	{
		final int particleBlock = 2 * particleId + 2;
		final int particleSize = blockStarts[particleBlock + 1] - blockStarts[particleBlock];
		final int particleStart = blockStarts[particleBlock];
		final Map<IntCoordinates, SparseMatrix> blocks = new HashMap<>();
		blocks.put(new IntCoordinates(0, 0),
		           zeroDerivativeMatrix.getBlocks()
		                               .get(new IntCoordinates(particleStart, particleStart))
		                               .add(twoDerivativeMatrix.getBlocks()
		                                                       .get(new IntCoordinates(particleStart,
		                                                                               particleStart))));
		blocks.put(new IntCoordinates(0, particleSize),
		           zeroDerivativeMatrix.getBlocks()
		                               .get(new IntCoordinates(particleStart, particleStart + particleSize)));
		blocks.put(new IntCoordinates(particleSize, 0),
		           oneDerivativeMatrix.getBlocks()
		                              .get(new IntCoordinates(particleStart + particleSize, particleStart)));
		return new BlockSparseMatrix(blocks, 2 * particleSize,
		                             2 * particleSize);
	}
	
	protected Matrix getAllBackgroundMatrices()
	{
		final Map<IntCoordinates, SparseMatrix> blocks = new HashMap<>();
		blocks.put(new IntCoordinates(0, 0),
		           zeroDerivativeMatrix.getBlocks()
		                               .get(new IntCoordinates(blockStarts[0], blockStarts[0]))
		                               .add(
			                               oneDerivativeMatrix.getBlocks()
			                                                  .get(new IntCoordinates(blockStarts[0],
			                                                                          blockStarts[0]))));
		blocks.put(new IntCoordinates(blockStarts[1] - blockStarts[0], 0),
		           zeroDerivativeMatrix.getBlocks()
		                               .get(new IntCoordinates(blockStarts[1], blockStarts[0])));
		blocks.put(new IntCoordinates(0, blockStarts[1] - blockStarts[0]),
		           zeroDerivativeMatrix.getBlocks()
		                               .get(new IntCoordinates(blockStarts[0], blockStarts[1])));
		return new BlockSparseMatrix(blocks, blockStarts[2] - blockStarts[0],
		                             blockStarts[2] - blockStarts[0]);
	}
	
	protected Matrix getAllBackgroundTransferMatrices(final int particleId)
	{
		final int particleBlock = 2 * particleId + 2;
		final int particleSize = blockStarts[particleBlock + 1] - blockStarts[particleBlock];
		final int particleStart = blockStarts[particleBlock];
		final Map<IntCoordinates, SparseMatrix> blocks = new HashMap<>();
		System.out.println(zeroDerivativeMatrix.getBlocks()
		                                       .keySet());
		System.out.println("aaaaa " + blockStarts[0] + " " + blockStarts[2 * particleId + 3]);
		blocks.put(new IntCoordinates(0, particleSize), evaluateBackgroundLagrangeIntegrals(particleId));
//		           zeroDerivativeMatrix.getBlocks()
//		                               .get(new IntCoordinates(blockStarts[0],
//		                                                       blockStarts[2 * particleId + 3])));
		return new BlockSparseMatrix(blocks, blockStarts[2] - blockStarts[0],
		                             2 * particleSize);
	}
	
	protected Matrix getAllBackgroundTransferTransposeMatrices(final int particleId)
	{
		final int particleBlock = 2 * particleId + 2;
		final int particleSize = blockStarts[particleBlock + 1] - blockStarts[particleBlock];
		final int particleStart = blockStarts[particleBlock];
//		return zeroDerivativeMatrix.slice(new IntCoordinates(0, particleStart),
//		                                  new IntCoordinates(blockStarts[2],
//		                                                     particleStart + 2 * particleSize))
//		                           .transpose();
		final Map<IntCoordinates, SparseMatrix> blocks = new HashMap<>();
		blocks.put(new IntCoordinates(particleSize, 0),
		           evaluateBackgroundLagrangeIntegrals(particleId).transpose());
//		           zeroDerivativeMatrix.getBlocks()
//		                               .get(new IntCoordinates(particleStart + particleSize, blockStarts[0])));
		return new BlockSparseMatrix(blocks, 2 * particleSize,
		                             blockStarts[2] - blockStarts[0]);
	}
	
	abstract protected Function<CoordinateVector, CoordinateVector> velocityBoundaryValues();
	
	protected Predicate<TPFace> getDirichletBoundary()
	{
		return f -> true;
	}
	
	@Override
	protected Function2<Matrix, MutableVector, Matrix> boundaryApplier()
	{
		return (mat, vec) ->
		{
			final MutableMatrix matrix;
			if (mat instanceof MutableMatrix)
				matrix = (MutableMatrix) mat;
			else if (mat.isSparse())
				matrix = new SparseMatrix(mat);
			else
				matrix = new DenseMatrix(mat);
			final VectorFunction velocityBoundaryFunction
				= VectorFunction.fromLambda(velocityBoundaryValues(),
				                            2, 2);
			backgroundSpace.writeBoundaryValuesTo(new ComposedMixedFunction(velocityBoundaryFunction),
			                                      getDirichletBoundary(),
			                                      (f, function) -> function.hasVelocityFunction(),
			                                      matrix,
			                                      vec);
			final Optional<ST> firstPressure =
				backgroundSpace.getShapeFunctions()
				               .values()
				               .stream()
				               .filter(ComposeMixedShapeFunction::hasPressureFunction)
				               .findFirst();
			firstPressure.ifPresent(st -> backgroundSpace.overWriteValue(st.getGlobalIndex(),
			                                                             0,
			                                                             matrix,
			                                                             vec));
			return matrix;
		};
	}
}
