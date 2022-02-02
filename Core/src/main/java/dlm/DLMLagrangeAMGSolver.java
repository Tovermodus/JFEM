package dlm;

import basic.PlotWindow;
import basic.VectorFunctionOnCells;
import io.vavr.Tuple2;
import it.unimi.dsi.fastutil.ints.Int2DoubleMap;
import linalg.*;
import mixed.*;
import multigrid.AMGPreconditionerSpace;
import multigrid.Smoother;
import org.jetbrains.annotations.NotNull;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

import java.util.ArrayList;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class DLMLagrangeAMGSolver
	extends DLMSolver
{
	public double smootherOmega = 1;
	public int smootherRuns = 7;
	
	PreconditionedIterativeImplicitSchur schur;
	MultiGridFluid fluid;
	final List<Particle> particles;
	
	public DLMLagrangeAMGSolver(final MultiGridFluid f, final List<Particle> particles)
	{
		fluid = f;
		this.particles = particles;
	}
	
	@NotNull
	private AMGPreconditionerSpace<TaylorHoodSpace, TPCell, TPFace, QkQkFunction, MixedValue, MixedGradient,
		MixedHessian> create_space(
		final List<ParticleIterate> particleStates, final double dt, final double t,
		final FluidIterate iterate)
	{
		return new AMGPreconditionerSpace<>(fluid.refinements, fluid.polynomialDegree)
		{
			@Override
			public List<TaylorHoodSpace> createSpaces(final int refinements)
			{
				final List<TaylorHoodSpace> ret = new ArrayList<>();
				for (int i = 0; i < refinements + 1; i++)
				{
					final TaylorHoodSpace subSpace = new TaylorHoodSpace(fluid.startCoordinates,
					                                                     fluid.endCoordinates,
					                                                     fluid.coarsestCells.mul((int) Math.pow(
						                                                     2,
						                                                     i)));
					subSpace.assembleCells();
					subSpace.assembleFunctions(fluid.polynomialDegree);
					ret.add(subSpace);
				}
				return ret;
			}
			
			@Override
			public Tuple2<SparseMatrix, DenseVector> createSystem(final TaylorHoodSpace space)
			{
				final int n = space.getShapeFunctionMap()
				                   .size();
				final VectorFunctionOnCells<TPCell, TPFace> velocity;
				final Vector restrictedIterate;
				if (iterate != null)
				{
					restrictedIterate = restrictToSize(n, iterate.current);
				} else
				{
					restrictedIterate = new DenseVector(n);
				}
				velocity = new MixedTPFESpaceFunction<>(space.getShapeFunctionMap(),
				                                        restrictedIterate).getVelocityFunction();
				final FluidSystem fs = fluid.getFluidSystemForSpace(space, velocity, 0,
				                                                    restrictedIterate);
				final var blockRhs = Fluid.getBlockRhs(fs, dt);
				final SparseMatrix s = new SparseMatrix(blockRhs._1);
				final List<Matrix> schurContributions =
					IntStream.range(0, particles.size())
					         .parallel()
					         .mapToObj(i -> new Tuple2<>(i,
					                                     particles.get(i)
					                                              .buildLagrangeBackgroundMatrix(
						                                              space,
						                                              particleStates.get(i))))
					         .map(e ->
					              {
						              final SparseMatrix block =
							              new SparseMatrix(s.getRows(),
							                               2 * e._2.getCols());
						              block.addSmallMatrixInPlaceAt(e._2,
						                                            0,
						                                            e._2.getCols());
						              return new Tuple2<>(e._1, block);
					              })
					         .map(e ->
						              e._2.mmMul(new SparseMatrix(schur.getDiagonalInverse(
							               e._1)))
						                  .mtMul(e._2)
						                  .mul(-1))
					         .collect(
						         Collectors.toList());
				for (final Matrix m : schurContributions)
					s.addInPlace(m);
				
				final Int2DoubleMap nodeValues = fluid.getDirichletNodeValuesForSpace(space, t);
				nodeValues.forEach((node, val) ->
				                   {
					                   s.deleteColumn(node);
					                   s.deleteRow(node);
					                   s.set(1, node, node);
				                   });
				return new Tuple2<>(s, new DenseVector(n));
			}
			
			@Override
			public List<Smoother> createSmoothers()
			{
				
				final List<Smoother> ret = new ArrayList<>();
				for (int i = 1; i < fluid.refinements + 1; i++)
				{
					ret.add(new BSSmoother2(smootherRuns,
					                        smootherOmega,
					                        spaces.get(i)
					                              .getVelocitySize()));
					final Int2DoubleMap nodeValues =
						fluid.getDirichletNodeValuesForSpace(spaces.get(i), t);
					final int finalI = i;
					nodeValues.forEach((node, val) ->
					                   {
						                   systems.get(finalI)
						                          .deleteColumn(node);
						                   systems.get(finalI)
						                          .deleteRow(node);
						                   systems.get(finalI)
						                          .set(1, node, node);
					                   });
				}
				final DenseVector v1 = new DenseVector(systems.get(0)
				                                              .getVectorSize());
				final DenseVector v2 = new DenseVector(systems.get(0)
				                                              .getVectorSize());
				int j = 6;
				for (final IntCoordinates c : v1.getShape()
				                                .range())
				{
					v1.set(j++, c);
					v2.set(1, c);
				}
				final SparseMatrix P = prolongationMatrices.get(0);
				final SparseMatrix PAP = new SparseMatrix(P.tmMul((SparseMatrix) systems.get(1))
				                                           .mmMul(P));
				final SparseMatrix s0 = new SparseMatrix(systems.get(0));
				final Int2DoubleMap nodeValues = fluid.getDirichletNodeValuesForSpace(spaces.get(0), t);
				nodeValues.forEach((node, val) ->
				                   {
					                   PAP.deleteColumn(node);
					                   PAP.deleteRow(node);
					                   PAP.set(1, node, node);
					                   s0.deleteColumn(node);
					                   s0.deleteRow(node);
					                   s0.set(1, node, node);
				                   });
				applyZeroBoundaryConditions(spaces.get(0), v1);
				applyZeroBoundaryConditions(spaces.get(0), v2);
				System.out.println(Math.abs(PAP.mvMul(v1)
				                               .inner(v2) - s0
					.mvMul(v1)
					.inner(v2)) + "   PAPPPP");
				
				PlotWindow.addPlot(new MatrixPlot(PAP,
				                                  "pap"));
				PlotWindow.addPlot(new MatrixPlot(((SparseMatrix) s0),
				                                  "s0"));
				PlotWindow.addPlot(new MatrixPlot((PAP).sub((SparseMatrix) s0),
				                                  "diff"));
				
				verbose = true;
				return ret;
			}
			
			@Override
			public void applyZeroBoundaryConditions(final TaylorHoodSpace space, final MutableVector vector)
			{
				final Int2DoubleMap nodeValues = fluid.getDirichletNodeValuesForSpace(space, t);
				nodeValues.forEach((node, val) ->
				                   {
					                   vector.set(0, node);
				                   });
			}
		};
	}
	
	@Override
	protected Vector solve(final BlockSparseMatrix systemMatrix,
	                       final DenseVector rhs,
	                       final FluidIterate fluidState,
	                       final List<ParticleIterate> particleStates,
	                       final FluidSystem fluidSystem,
	                       final List<ParticleSystem> particleSystems, final double dt, final double t)
	{
		
		if (schur == null)
			schur = new PreconditionedIterativeImplicitSchur(systemMatrix,
			                                                 null);
		else
			schur.resetOffDiagonals(systemMatrix);
		final AMGPreconditionerSpace<?, ?, ?, ?, ?, ?, ?> mg = create_space(particleStates, dt, t, fluidState);
		final DenseVector d = new DenseVector(mg.getVectorSize());
		for (int i = 0; i < d.getLength(); i++)
			d.set(Math.random(), i);
//		final ExplicitSchurSolver sss = new DirectSchur(systemMatrix);
//		PlotWindow.addPlot(new MatrixPlot(sss.getSchurComplement(), "true"));
//		PlotWindow.addPlot(new MatrixPlot(new DenseMatrix((SparseMatrix) mg.finest_system), "mg"));
//		PlotWindow.addPlot(new MatrixPlot(((SparseMatrix) mg.finest_system).sub(new SparseMatrix(sss.getSchurComplement())),
//		                                  "diff"));
//
		System.out.println(schur.schurMvMul()
		                        .mvMul(d)
		                        .sub(mg.finest_system.mvMul(d))
		                        .absMaxElement() + "  diffffaskdjfh");
		schur.preconditioner = mg;
		return schur.mvMul(rhs);
	}
}
