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
		systems = new ArrayList<>();
		int k = 0;
		for (final CSpace space : spaces)
		{
			space.assembleCells();
			space.assembleFunctions(polynomialDegree);
			final var mat_rhs = createSystem(space);
			systems.add(mat_rhs._1);
			finest_rhs = mat_rhs._2;
			finest_system = mat_rhs._1;
			System.out.println("functions, system " + k++);
		}
		smoothers = createSmoothers();
		System.out.println("smoothers");
		prolongationOperators = new ArrayList<>();
		prolongationMatrices = new ArrayList<>();
		for (int i = 0; i < spaces.size() - 1; i++)
			prolongationOperators.add(createProlongationMatrix(spaces.get(i), spaces.get(i + 1)));
		System.out.println("prolongation");
	}
	
	SparseMvMul createProlongationMatrix(final CSpace coarse, final CSpace fine)
	{
		final SparseMatrix prolongationMatrix = new SparseMatrix(fine.getShapeFunctions()
		                                                             .size(),
		                                                         coarse.getShapeFunctions()
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
						                                    coarse.getCellSupportMapping()
						                                          .get(
							                                          coarseCell);
					                                    final Collection<ST> fineFunctions =
						                                    fine.getCellSupportMapping()
						                                        .get(fineCell);
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
		//System.out.println(guess.getLength());
		if (level == 0)
		{
			if (systems.get(0) instanceof Matrix)
				return new DenseMatrix((Matrix) systems.get(0)).solve(rhs);
			final Vector solution = new IterativeSolver(true).solveGMRES(systems.get(0), rhs, 1e-9);
			return solution;
		}
		guess = smoothers.get(level - 1)
		                 .smooth(systems.get(level), rhs, guess);
		final Vector restrictedDefect =
			prolongationOperators.get(level - 1)
			                     .tvMul(rhs.sub(systems.get(level)
			                                           .mvMul(guess)));
		final Vector correction = mgStep(level - 1,
		                                 new DenseVector(restrictedDefect.getLength()),
		                                 restrictedDefect);
		guess = guess.add(prolongationOperators.get(level - 1)
		                                       .mvMul(correction));
		guess = smoothers.get(level - 1)
		                 .smooth(systems.get(level), rhs, guess);
		return guess;
	}
	
	@Override
	public Vector vCycle(final Vector initialIterate, final Vector rhs)
	{
		return mgStep(spaces.size() - 1, initialIterate, rhs);
	}
	
	@Override
	public VectorMultiplyable getFinestSystem()
	{
		return finest_system;
	}
}
