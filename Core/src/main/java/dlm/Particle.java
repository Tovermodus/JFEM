package dlm;

import basic.*;
import com.google.common.base.Stopwatch;
import distorted.*;
import distorted.geometry.DistortedCell;
import io.vavr.Tuple2;
import it.unimi.dsi.fastutil.ints.Int2DoubleMap;
import linalg.*;
import mixed.ComposeMixedShapeFunction;
import mixed.MixedPlot2DTime;
import mixed.TaylorHoodSpace;
import org.jetbrains.annotations.NotNull;
import tensorproduct.ContinuousTPVectorFunction;

import java.util.List;
import java.util.function.Function;

public interface Particle
{
	
	DistortedGridSpace<DistortedVectorShapeFunction, CoordinateVector, CoordinateMatrix, CoordinateTensor> getSpace();
	
	List<CellIntegral<DistortedCell, DistortedVectorShapeFunction>> getIntegrals();
	
	List<CellIntegral<DistortedCell, DistortedVectorShapeFunction>> getLagrangeIntegrals();
	
	List<CellIntegral<DistortedCell, DistortedVectorShapeFunction>> getMassIntegrals();
	
	List<RightHandSideIntegral<DistortedCell, DistortedVectorShapeFunction>> getForceIntegrals(double t);
	
	List<RightHandSideIntegral<DistortedCell, DistortedVectorShapeFunction>> getBackgroundLagrangeIntegrals(
		DistortedVectorFunctionOnCells backgroundFunctionAtX);
	
	List<CellIntegral<DistortedCell, DistortedVectorShapeFunction>>
	getSemiImplicitIntegrals(final DistortedVectorFunctionOnCells displacement);
	
	default DistortedVectorFunctionOnCells getDisplacement(final ParticleIterate iterate)
	{
		return new DistortedVectorFESpaceFunction(getSpace().getShapeFunctionMap(), iterate.current);
	}
	
	default DistortedVectorFunctionOnCells getPosition(final ParticleIterate iterate)
	{
		final DistortedVectorFunctionOnCells displacement = getDisplacement(iterate);
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
				           .add(displacement.valueOnReferenceCell(pos, cell));
			}
			
			@Override
			public CoordinateMatrix gradientOnReferenceCell(final CoordinateVector pos,
			                                                final DistortedCell cell)
			{
				return CoordinateDenseMatrix.identity(2)
				                            .add(displacement.gradientOnReferenceCell(pos, cell));
			}
			
			@Override
			public int getDomainDimension()
			{
				return 2;
			}
			
			@Override
			public CoordinateVector value(final CoordinateVector pos)
			{
				return pos.add(displacement.value(pos));
			}
			
			@Override
			public CoordinateMatrix gradient(final CoordinateVector pos)
			{
				return CoordinateDenseMatrix.identity(2)
				                            .add(displacement.gradient(pos));
			}
		};
	}
	
	Function<CoordinateVector, CoordinateVector> getInitialVelocity();
	
	default ParticleIterate buildInitialIterate(final double dt)
	{
		final DenseVector initialDisplacement = new DenseVector(getSystemSize());
		final DenseVector initialLagrange = new DenseVector(getLagrangeSize());
		final DenseVector lastLagrange = new DenseVector(getLagrangeSize());
		final DenseVector initialVelocity = new DenseVector(getSystemSize());
		final VectorFunction initialVelo = VectorFunction.fromLambda(getInitialVelocity(), 2, 2);
		getSpace().getShapeFunctionMap()
		          .values()
		          .forEach(function ->
		                   {
			                   initialVelocity.set(function.getNodeFunctional()
			                                               .evaluate(initialVelo),
			                                       function.getGlobalIndex());
		                   });
		final DenseVector lastDisplacement = initialVelocity.mul(-dt);
		return new ParticleIterate(initialDisplacement, initialLagrange, lastDisplacement, lastLagrange);
	}
	
	default ParticleSystem buildSystem(final Fluid f,
	                                   final double t, final ParticleIterate iterate)
	{
		final Matrix massMatrix = buildMassMatrix();
		final Matrix elasticityMatrix = buildElasticityMatrix();
		final Matrix semiImplicitMatrix = buildSemiImplicitMatrix(iterate);
		final Matrix lagrangeMatrix = buildLagrangeMatrix();
		final Matrix lagrangeBackgroundMatrix = buildLagrangeBackgroundMatrix(f.getSpace(), iterate);
		final Vector forceRhs = buildForceRhs(t);
		final Vector accelerationRhs = buildAccelerationRhs(iterate, massMatrix);
		final Vector lagrangeRhs = buildLagrangeRhs(iterate, lagrangeMatrix);
		return new ParticleSystem(massMatrix, elasticityMatrix, semiImplicitMatrix,
		                          lagrangeMatrix,
		                          lagrangeBackgroundMatrix, forceRhs, accelerationRhs, lagrangeRhs);
	}
	
	default Tuple2<SparseMatrix, DenseVector> getBlockRhs(final ParticleSystem ps, final double dt)
	{
		final SparseMatrix displacementMatrix =
			new SparseMatrix(
				ps.massMatrix.mul(1. / (dt))
				             .add(ps.elasticityMatrix.mul(dt))
				             .add(ps.semiImplicitMatrix.mul(dt)));
		final DenseVector displacementRhs =
			new DenseVector(ps.forceRhs.mul(1)
			                           .add(ps.accelerationRhs.mul(1. / (dt * dt))));
		
		final SparseMatrix s = new SparseMatrix(getSystemSize() + getLagrangeSize(),
		                                        getSystemSize() + getLagrangeSize());
		s.addSmallMatrixInPlaceAt(displacementMatrix, 0, 0);
		s.addSmallMatrixInPlaceAt(ps.lagrangeMatrix.mul(1), 0, getSystemSize());
		s.addSmallMatrixInPlaceAt(ps.lagrangeMatrix.transpose()
		                                           .mul(1.), getSystemSize(), 0);
		final DenseVector d = new DenseVector(getSystemSize() + getLagrangeSize());
		d.addSmallVectorAt(displacementRhs, 0);
		d.addSmallVectorAt(ps.lagrangeRhs.mul(1. / dt), getSystemSize());
		return new Tuple2<>(s, d);
	}
	
	default SparseMatrix getLagrangeBackgroundBlock(final ParticleSystem ps, final Fluid f,
	                                                final double dt)
	{
		final SparseMatrix s = new SparseMatrix(f.getSystemSize(), getSystemSize() + getLagrangeSize());
		s.addSmallMatrixInPlaceAt(ps.lagrangeBackgroundMatrix.mul(1), 0, getSystemSize());
		return s;
	}
	
	default SparseMatrix getLagrangeBackgroundBlockTranspose(final ParticleSystem ps,
	                                                         final Fluid f,
	                                                         final double dt)
	{
		return getLagrangeBackgroundBlock(ps, f, dt).transpose();
	}
	
	private static Vector buildLagrangeRhs(final ParticleIterate iterate,
	                                       final Matrix lagrangeMatrix)
	{
		return lagrangeMatrix.tvMul(iterate.current);
	}
	
	private Vector buildForceRhs(final double t)
	{
		final DenseVector rhs = new DenseVector(getSystemSize());
		getSpace().writeCellIntegralsToRhs(getForceIntegrals(t), rhs);
		return rhs;
	}
	
	private static Vector buildAccelerationRhs(final ParticleIterate iterate,
	                                           final Matrix massMatrix)
	{
		return massMatrix.mvMul(iterate.current
			                        .mul(2.)
			                        .sub(iterate.last));
	}
	
	default Matrix buildLagrangeBackgroundMatrix(final TaylorHoodSpace fluidSpace, final ParticleIterate iterate)
	{
		
		final CountMetric cm = MetricWindow.getInstance()
		                                   .setMetric("lagBackGround",
		                                              new CountMetric(getSpace().getCells()
		                                                                        .size()));
		final SparseMatrix lagrangeBackgroundMatrix = new SparseMatrix(fluidSpace.getVelocitySize(),
		                                                               getLagrangeSize());
		final DistortedVectorFunctionOnCells X = getPosition(iterate);
		final double particleCellDiam = getSpace().getMaxDiam();
		final double fluidCellDiam = fluidSpace
			.getMaxDiam();
		System.out.println("lagback");
		final Stopwatch s = Stopwatch.createStarted();
		
		getSpace().forEachCell(cell ->
		                       {
			                       cm.increment();
			                       fluidSpace
				                       .getShapeFunctionMap()
				                       .values()
				                       .stream()
				                       .parallel()
				                       .filter(ComposeMixedShapeFunction::hasVelocityFunction)
				                       .map(ComposeMixedShapeFunction::getVelocityFunction)
				                       .filter(backGroundFunction -> hasSupportOverlap(X,
				                                                                       particleCellDiam,
				                                                                       fluidCellDiam,
				                                                                       cell,
				                                                                       backGroundFunction))
				                       .forEach(backgroundFunction ->
					                                writeBackgroundLagrangeIntegralsOnCellToMatrix(
						                                lagrangeBackgroundMatrix,
						                                cell,
						                                backgroundFunction,
						                                X));
		                       });
		System.out.println("lagbackDone" + s.elapsed());
		return lagrangeBackgroundMatrix;
	}
	
	private static boolean hasSupportOverlap(final DistortedVectorFunctionOnCells X,
	                                         final double particleCellDiam,
	                                         final double fluidCellDiam,
	                                         final DistortedCell cell,
	                                         final ContinuousTPVectorFunction backGroundFunction)
	{
		return backGroundFunction.getNodeFunctionalPoint()
		                         .sub(X.valueOnReferenceCell(cell.referenceCell.center(), cell))
		                         .euclidianNorm() < 2 * (particleCellDiam + fluidCellDiam);
	}
	
	private void writeBackgroundLagrangeIntegralsOnCellToMatrix(final SparseMatrix lagrangeBackGroundMatrix,
	                                                            final DistortedCell cell,
	                                                            final ContinuousTPVectorFunction backgroundFunction,
	                                                            final DistortedVectorFunctionOnCells X)
	{
		final DistortedVectorFunctionOnCells backgroundAtX =
			DistortedVectorFunctionOnCells.concatenate(backgroundFunction, X);
		final int rowIndex
			= backgroundFunction.getGlobalIndex();
		final var backGroundIntegrals = getBackgroundLagrangeIntegrals(backgroundAtX);
		getSpace().getShapeFunctionsWithSupportOnCell(cell)
		          .forEach(particleFunction ->
		                   {
			                   final int colIndex = particleFunction.getGlobalIndex();
			                   double value = 0;
			                   for (final var integral : backGroundIntegrals)
				                   value += integral.evaluateRightHandSideIntegral(cell,
				                                                                   particleFunction);
			                   if (value != 0)
				                   lagrangeBackGroundMatrix.add(value, rowIndex, colIndex);
		                   });
	}
	
	private Matrix buildLagrangeMatrix()
	{
		final SparseMatrix lagrangeMatrix = new SparseMatrix(getSystemSize(), getLagrangeSize());
		getSpace().writeCellIntegralsToMatrix(getLagrangeIntegrals(), lagrangeMatrix);
		return lagrangeMatrix.mul(-1);
	}
	
	private Matrix buildSemiImplicitMatrix(final ParticleIterate iterate)
	{
		final SparseMatrix semiImplicitMatrix = new SparseMatrix(getSystemSize(), getSystemSize());
		getSpace().writeCellIntegralsToMatrix(getSemiImplicitIntegrals(getDisplacement(iterate)),
		                                      semiImplicitMatrix);
		return semiImplicitMatrix;
	}
	
	@NotNull
	private Matrix buildMassMatrix()
	{
		final SparseMatrix massMatrix;
		massMatrix = new SparseMatrix(getSystemSize(), getSystemSize());
		getSpace().writeCellIntegralsToMatrix(getMassIntegrals(), massMatrix);
		return massMatrix;
	}
	
	@NotNull
	private Matrix buildElasticityMatrix()
	{
		final SparseMatrix elasticityMatrix;
		elasticityMatrix = new SparseMatrix(getSystemSize(), getSystemSize());
		getSpace().writeCellIntegralsToMatrix(getIntegrals(), elasticityMatrix);
		return elasticityMatrix;
	}
	
	default int getSystemSize()
	{
		return getSpace().getShapeFunctionMap()
		                 .size();
	}
	
	default int getLagrangeSize()
	{
		return getSpace().getShapeFunctionMap()
		                 .size();
	}
	
	default Overlay generateOverlay(final DenseMatrix particleHistory, final MixedPlot2DTime backgroundPlot)
	{
		return new DistortedOverlay(backgroundPlot,
		                            getSpace(),
		                            particleHistory,
		                            5);
	}
	
	Int2DoubleMap getDirichletNodeValues(double t);
}
