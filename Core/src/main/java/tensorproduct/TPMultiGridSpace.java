package tensorproduct;

import basic.ShapeFunction;
import com.google.common.collect.TreeMultimap;
import io.vavr.Tuple2;
import linalg.*;
import multigrid.Smoother;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;

public abstract class TPMultiGridSpace<CSpace extends CartesianGridSpace<ST, valueT, gradientT, hessianT>,
	ST extends ShapeFunction<TPCell, TPFace, valueT, gradientT, hessianT> & Comparable<ST>, valueT, gradientT,
	hessianT>
{
	protected final CoordinateVector startCoordinates;
	protected final CoordinateVector endCoordinates;
	protected final IntCoordinates coarseCellsPerDimension;
	public List<CSpace> spaces;
	public List<Smoother> smoothers;
	public List<SparseMvMul> prolongationOperators;
	public List<SparseMatrix> prolongationMatrices;
	public List<VectorMultiplyable> systems;
	public Vector finest_rhs;
	public VectorMultiplyable finest_system;
	
	public TPMultiGridSpace(final CoordinateVector startCoordinates, final CoordinateVector endCoordinates,
	                        final IntCoordinates coarseCellsPerDimension, final int refinements,
	                        final int polynomialDegree)
	{
		this.startCoordinates = startCoordinates;
		this.endCoordinates = endCoordinates;
		this.coarseCellsPerDimension = coarseCellsPerDimension;
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
					                                    final Set<ST> coarseFunctions =
						                                    coarse.supportOnCell.get(
							                                    coarseCell);
					                                    final Set<ST> fineFunctions =
						                                    fine.supportOnCell.get(
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
	
	public Vector mgStep(final int level, Vector guess, final Vector rhs)
	{
		//System.out.println("mgstep" + level);
		//System.out.println(guess.getLength());
		if (level == 0)
		{
			final Vector solution = new IterativeSolver(false).solveGMRES(systems.get(0), rhs, 1e-9);
			return solution;
		}
		for (int i = 0; i < 2; i++)
			guess = smoothers.get(level - 1)
			                 .smooth(systems.get(level), rhs, guess);
		Vector restrictedDefect =
			prolongationOperators.get(level - 1)
			                     .tvMul(rhs.sub(systems.get(level)
			                                           .mvMul(guess)));
		final Vector correction = mgStep(level - 1,
		                                 new DenseVector(restrictedDefect.getLength()),
		                                 restrictedDefect);
		guess = guess.add(prolongationOperators.get(level - 1)
		                                       .mvMul(correction));
		restrictedDefect =
			prolongationOperators.get(level - 1)
			                     .tvMul(rhs.sub(systems.get(level)
			                                           .mvMul(guess)));
		
		for (int i = 0; i < 2; i++)
			guess = smoothers.get(level - 1)
			                 .smooth(systems.get(level), rhs, guess);
		//System.out.println("mgstep done" + level);
		return guess;
	}
	
	public Vector vCycle(final Vector initialIterate, final Vector rhs)
	{
		return mgStep(spaces.size() - 1, initialIterate, rhs);
	}
}
