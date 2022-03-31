package multigrid;

import basic.*;
import com.google.common.collect.TreeMultimap;
import io.vavr.Tuple2;
import linalg.*;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

public abstract class MGSpace<CSpace extends FESpace<CT, FT, ST> & Assembleable, CT extends Cell<CT, FT>,
	FT extends Face<CT, FT>,
	ST extends ShapeFunction<CT, FT, valueT, gradientT, hessianT> & Comparable<ST>, valueT, gradientT, hessianT>
	implements MGInterface
{
	public List<CSpace> spaces;
	public List<Smoother> smoothers;
	public List<SparseMvMul> prolongationOperators;
	public List<SparseMatrix> prolongationMatrices;
	public List<VectorMultiplyable> systems;
	public Vector finest_rhs;
	public VectorMultiplyable finest_system;
	
	public MGSpace(final int refinements,
	               final int polynomialDegree)
	{
		spaces = createSpaces(refinements);
		System.out.println("spaces");
		prolongationOperators = new ArrayList<>();
		prolongationMatrices = new ArrayList<>();
		for (final CSpace space : spaces)
		{
			space.assembleCells();
			space.assembleFunctions(polynomialDegree);
		}
		for (int i = 0; i < spaces.size() - 1; i++)
			prolongationOperators.add(createProlongationMatrix(spaces.get(i), spaces.get(i + 1)));
		System.out.println("prolongation");
		systems = new ArrayList<>();
		int k = 0;
		for (final CSpace space : spaces)
		{
			final var mat_rhs = createSystem(space);
			systems.add(mat_rhs._1);
			finest_rhs = mat_rhs._2;
			finest_system = mat_rhs._1;
			System.out.println("functions, system " + k++);
		}
		smoothers = createSmoothers();
		System.out.println("smoothers");
	}
	
	public Vector restrictToSize(final int newSize, Vector finest)
	{
		final int startLength = finest.getLength();
		if (spaces.stream()
		          .allMatch(s -> s.getShapeFunctionMap()
		                          .size() != startLength))
			throw new IllegalArgumentException("newsize not valid");
		int l = prolongationOperators.size() - 1;
		while (finest.getLength() != newSize)
		{
			finest = new DenseVector(prolongationOperators.get(l)
			                                              .tvMul(finest));
			l--;
		}
		return finest;
	}
	
	SparseMvMul createProlongationMatrix(final CSpace coarse, final CSpace fine)
	{
		final SparseMatrix prolongationMatrix = new SparseMatrix(fine.getShapeFunctionMap()
		                                                             .size(),
		                                                         coarse.getShapeFunctionMap()
		                                                               .size());
		final TreeMultimap<ST, ST> refinedFunctions = getRefinedFunctions(coarse, fine);
		refinedFunctions.forEach((coarseFunction, fineFunction) ->
		                         {
			                         prolongationMatrix.add(fineFunction.getNodeFunctional()
			                                                            .evaluate(coarseFunction),
			                                                fineFunction.getGlobalIndex(),
			                                                coarseFunction.getGlobalIndex());
		                         });
		prolongationMatrices.add(prolongationMatrix);
		return new SparseMvMul(prolongationMatrix);
	}
	
	private TreeMultimap<ST, ST> getRefinedFunctions(final CSpace coarse, final CSpace fine)
	{
		final TreeMultimap<ST, ST> ret = TreeMultimap.create();
		coarse.forEachCell(coarseCell ->
			                   fine.forEachCell(fineCell ->
			                                    {
				                                    if (coarseCell.isInCell(fineCell.center()))
				                                    {
					                                    final Collection<ST> coarseFunctions =
						                                    coarse.getShapeFunctionsWithSupportOnCell(
							                                    coarseCell);
					                                    final Collection<ST> fineFunctions =
						                                    fine.getShapeFunctionsWithSupportOnCell(
							                                    fineCell);
					                                    synchronized (this)
					                                    {
						                                    for (final ST function : coarseFunctions)
						                                    {
							                                    ret.putAll(
								                                    function,
								                                    fineFunctions);
						                                    }
					                                    }
				                                    }
			                                    }));
		return ret;
	}
	
	public abstract List<CSpace> createSpaces(int refinements);
	
	public abstract Tuple2<VectorMultiplyable, DenseVector> createSystem(CSpace space);
	
	public abstract List<Smoother> createSmoothers();
	
	@Override
	public Vector mgStep(final int level, Vector guess, final Vector rhs)
	{
		if (level == 0)
		{
			System.out.println("    level 0 " + systems.get(0)
			                                           .getVectorSize());
			if (systems.get(0) instanceof Matrix)
			{
				final DenseVector solve = new DenseMatrix((Matrix) systems.get(0)).solve(rhs);
				System.out.println("    solved");
				return solve;
			}
			final Vector solution = new IterativeSolver(true).solveGMRES(systems.get(0), rhs, 1e-9);
			return solution;
		}
		System.out.println(" MG level " + level + " residual " + rhs.sub(systems.get(level)
		                                                                        .mvMul(guess))
		                                                            .absMaxElement());
		guess = smoothers.get(level - 1)
		                 .smooth(systems.get(level), rhs, guess);
		Vector restrictedDefect =
			prolongationOperators.get(level - 1)
			                     .tvMul(rhs.sub(systems.get(level)
			                                           .mvMul(guess)));
		System.out.println(" MG level " + level + " precorrection residual " + rhs.sub(systems.get(level)
		                                                                                      .mvMul(guess))
		                                                                          .absMaxElement());
		System.out.println(restrictedDefect.absMaxElement());
		final Vector correction = mgStep(level - 1,
		                                 new DenseVector(restrictedDefect.getLength()),
		                                 restrictedDefect);
		guess = guess.add(prolongationOperators.get(level - 1)
		                                       .mvMul(correction));
		restrictedDefect =
			prolongationOperators.get(level - 1)
			                     .tvMul(rhs.sub(systems.get(level)
			                                           .mvMul(guess)));
		System.out.println(restrictedDefect.absMaxElement());
		System.out.println(" MG level " + level + " correction residual " + rhs.sub(systems.get(level)
		                                                                                   .mvMul(guess))
		                                                                       .absMaxElement());
		guess = smoothers.get(level - 1)
		                 .smooth(systems.get(level), rhs, guess);
		System.out.println(" MG level " + level + " final residual " + rhs.sub(systems.get(level)
		                                                                              .mvMul(guess))
		                                                                  .absMaxElement());
		return guess;
	}
	
	@Override
	public Vector vCycle(final Vector initialIterate, final Vector rhs)
	{
//		try
//		{
//			Thread.sleep(100);
//		} catch (final InterruptedException e)
//		{
//			e.printStackTrace();
//		}
		return mgStep(spaces.size() - 1, initialIterate, rhs);
	}
	
	@Override
	public Vector fullVCycleSolver(Vector rhs)
	{
		final List<Vector> rhsides = new ArrayList<>();
		rhsides.add(rhs);
		for (int i = prolongationOperators.size() - 1; i >= 0; i--)
		{
			rhs = prolongationOperators.get(i)
			                           .tvMul(rhs);
			rhsides.add(rhs);
		}
		Vector iterate = new DenseVector(systems.get(0)
		                                        .getVectorSize());
		for (int i = 0; i < prolongationOperators.size(); i++)
		{
			iterate = mgStep(i, iterate, rhsides.get(spaces.size() - i - 1));
			iterate = prolongationOperators.get(i)
			                               .mvMul(iterate);
		}
		return mgStep(spaces.size() - 1, iterate, rhsides.get(0));
	}
	
	@Override
	public VectorMultiplyable getFinestSystem()
	{
		return finest_system;
	}
}
