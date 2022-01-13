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

public abstract class TPAlgebraicMultiGridSpace<CSpace extends CartesianGridSpace<ST, valueT, gradientT, hessianT>,
	ST extends ShapeFunction<TPCell, TPFace, valueT, gradientT, hessianT> & Comparable<ST>, valueT, gradientT,
	hessianT>
{
	protected final CoordinateVector startCoordinates;
	protected final CoordinateVector endCoordinates;
	protected final IntCoordinates coarseCellsPerDimension;
	public List<CSpace> spaces;
	public List<Smoother> smoothers;
	public List<SparseMvMul> prolongationOperator;
	public List<SparseMatrix> prolongationMatrices;
	public List<SparseMatrix> systems;
	public List<SparseMvMul> restrictionOperator;
	public List<SparseMatrix> restrictionMatrices;
	public Vector rhs;
	public SparseMatrix matrix;
	
	public TPAlgebraicMultiGridSpace(final CoordinateVector startCoordinates, final CoordinateVector endCoordinates,
	                                 final IntCoordinates coarseCellsPerDimension, final int refinements,
	                                 final int polynomialDegree)
	{
		this.startCoordinates = startCoordinates;
		this.endCoordinates = endCoordinates;
		this.coarseCellsPerDimension = coarseCellsPerDimension;
		spaces = createSpaces(refinements);
		System.out.println("spaces");
		for (final CSpace space : spaces)
		{
			space.assembleCells();
			space.assembleFunctions(polynomialDegree);
		}
		System.out.println("functions");
		final Tuple2<SparseMatrix, Vector> topLevelSystem =
			createFinestLevelSystem(spaces.get(spaces.size() - 1));
		System.out.println("finest system");
		matrix = topLevelSystem._1;
		rhs = topLevelSystem._2;
		prolongationMatrices = new ArrayList<>();
		restrictionMatrices = new ArrayList<>();
		prolongationOperator = new ArrayList<>();
		restrictionOperator = new ArrayList<>();
		for (int i = 0; i < spaces.size() - 1; i++)
			prolongationOperator.add(createProlongationMatrix(spaces.get(i), spaces.get(i + 1)));
		System.out.println("prolongation");
		for (int i = 0; i < spaces.size() - 1; i++)
			restrictionOperator.add(createRestrictionMatrix(spaces.get(i), spaces.get(i + 1)));
		System.out.println("restriction");
		System.out.println();
		systems = createSystems();
		System.out.println("systems");
		smoothers = createSmoothers();
		System.out.println("smoothers");
	}
	
	private List<SparseMatrix> createSystems()
	{
		final List<SparseMatrix> reversed = new ArrayList<>();
		reversed.add(matrix);
		for (int i = 1; i < spaces.size(); i++)
		{
			System.out.println("system " + i);
			final int systemIndex = spaces.size() - i - 1;
			final SparseMatrix coarser =
				new SparseMatrix(restrictionMatrices.get(systemIndex)
				                                    .mmMul(reversed.get(i - 1))
				                                    .mmMul(prolongationMatrices.get(systemIndex)));
			reversed.add(coarser);
		}
		final List<SparseMatrix> ret = new ArrayList<>();
		for (int i = spaces.size() - 1; i >= 0; i--)
			ret.add(reversed.get(i));
		return ret;
	}
	
	SparseMvMul createProlongationMatrix(final CSpace coarse, final CSpace fine)
	{
		final SparseMatrix prolongationMatrix = new SparseMatrix(fine.getShapeFunctions()
		                                                             .size(),
		                                                         coarse.getShapeFunctions()
		                                                               .size());
		final TreeMultimap<ST, ST> refinedFunctions = getRefinedFunctions(coarse, fine);
		refinedFunctions.forEach((coarseFunction, fineFunction) ->
			                         prolongationMatrix.add(fineFunction.getNodeFunctional()
			                                                            .evaluate(coarseFunction),
			                                                fineFunction.getGlobalIndex(),
			                                                coarseFunction.getGlobalIndex()));
		prolongationMatrices.add(prolongationMatrix);
		return new SparseMvMul(prolongationMatrix);
	}
	
	SparseMvMul createRestrictionMatrix(final CSpace coarse, final CSpace fine)
	{
		final SparseMatrix restrictionMatrix = new SparseMatrix(coarse.getShapeFunctions()
		                                                              .size(),
		                                                        fine.getShapeFunctions()
		                                                            .size());
		final TreeMultimap<ST, ST> refinedFunctions = getRefinedFunctions(coarse, fine);
		refinedFunctions.forEach((coarseFunction, fineFunction) ->
			                         restrictionMatrix.add(coarseFunction.getNodeFunctional()
			                                                             .evaluate(fineFunction),
			                                               coarseFunction.getGlobalIndex(),
			                                               fineFunction.getGlobalIndex()));
		restrictionMatrices.add(restrictionMatrix);
		return new SparseMvMul(restrictionMatrix);
	}
	
	private TreeMultimap<ST, ST> getRefinedFunctions(final CSpace coarse, final CSpace fine)
	{
		final TreeMultimap<ST, ST> ret = TreeMultimap.create();
		coarse.forEachCell(coarseCell ->
		                   {
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
			                                    });
		                   });
		return ret;
	}
	
	public abstract List<CSpace> createSpaces(int refinements);
	
	public abstract Tuple2<SparseMatrix, Vector> createFinestLevelSystem(CSpace space);
	
	public abstract List<Smoother> createSmoothers();
	
	public Vector mgStep(final int level, Vector guess, final Vector rhs)
	{
		//System.out.println("mgstep" + level);
		//System.out.println(guess.getLength());
		if (level == 0)
		{
			final Vector solution = new IterativeSolver().solveGMRES(systems.get(0), rhs, 1e-9);
			return solution;
		}
		guess = smoothers.get(level - 1)
		                 .smooth(systems.get(level), rhs, guess);
		final Vector restrictedDefect =
			restrictionOperator.get(level - 1)
			                   .mvMul(rhs.sub(systems.get(level)
			                                         .mvMul(guess)));
		final Vector correction = mgStep(level - 1,
		                                 new DenseVector(restrictedDefect.getLength()),
		                                 restrictedDefect);
		guess = guess.add(prolongationOperator.get(level - 1)
		                                      .mvMul(correction));
		
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