package tensorproduct;

import basic.DoubleCompare;
import basic.ShapeFunction;
import com.google.common.collect.TreeMultimap;
import io.vavr.Tuple2;
import linalg.*;
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
	public List<VectorMultiplyable> smoothers;
	public List<SparseMvMul> prolongationMatrices;
	public List<SparseMatrix> prolongationMatrices2;
	public List<Tuple2<VectorMultiplyable, DenseVector>> systems;
	public List<SparseMvMul> restrictionMatrices;
	
	public TPMultiGridSpace(final CoordinateVector startCoordinates, final CoordinateVector endCoordinates,
	                        final IntCoordinates coarseCellsPerDimension, final int refinements,
	                        final int polynomialDegree)
	{
		this.startCoordinates = startCoordinates;
		this.endCoordinates = endCoordinates;
		this.coarseCellsPerDimension = coarseCellsPerDimension;
		spaces = createSpaces(refinements);
		for (final CSpace space : spaces)
		{
			space.assembleCells();
			space.assembleFunctions(polynomialDegree);
		}
		systems = createSystems();
		smoothers = createSmoothers();
		prolongationMatrices = new ArrayList<>();
		prolongationMatrices2 = new ArrayList<>();
		for (int i = 0; i < spaces.size() - 1; i++)
			prolongationMatrices.add(createProlongationMatrix(spaces.get(i), spaces.get(i + 1)));
		restrictionMatrices = new ArrayList<>();
		for (int i = 0; i < spaces.size() - 1; i++)
			restrictionMatrices.add(createRestrictionMatrix(spaces.get(i), spaces.get(i + 1)));
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
			                         if (!fine.fixedNodes.contains(fineFunction.getGlobalIndex())
				                         && !coarse.fixedNodes.contains(coarseFunction.getGlobalIndex()))
				                         prolongationMatrix.add(fineFunction.getNodeFunctional()
				                                                            .evaluate(coarseFunction),
				                                                fineFunction.getGlobalIndex(),
				                                                coarseFunction.getGlobalIndex());
			                         else if (DoubleCompare.almostEqual(fineFunction.getNodeFunctional()
			                                                                        .evaluate(coarseFunction),
			                                                            1))
				                         prolongationMatrix.add(1, fineFunction.getGlobalIndex(),
				                                                coarseFunction.getGlobalIndex());
		                         });
		prolongationMatrices2.add(prolongationMatrix);
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
		                         {
			                         if (!fine.fixedNodes.contains(fineFunction.getGlobalIndex())
				                         && !coarse.fixedNodes.contains(coarseFunction.getGlobalIndex()))
				                         restrictionMatrix.add(coarseFunction.getNodeFunctional()
				                                                             .evaluate(fineFunction),
				                                               coarseFunction.getGlobalIndex(),
				                                               fineFunction.getGlobalIndex());
			                         else if (DoubleCompare.almostEqual(fineFunction.getNodeFunctional()
			                                                                        .evaluate(coarseFunction),
			                                                            1))
				                         restrictionMatrix.add(1, coarseFunction.getGlobalIndex(),
				                                               fineFunction.getGlobalIndex());
		                         });
		//PlotWindow.addPlot(new MatrixPlot(prolongationMatrix));
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
	
	public abstract List<Tuple2<VectorMultiplyable, DenseVector>> createSystems();
	
	public abstract List<VectorMultiplyable> createSmoothers();
	
	public Vector mgStep(final int level, Vector guess, final Vector rhs)
	{
		//System.out.println("mgstep" + level);
		//System.out.println(guess.getLength());
		if (level == 0)
		{
			final Vector solution = new IterativeSolver().solveGMRES(systems.get(0)._1, rhs, 1e-9);
			return solution;
		}
		for (int i = 0; i < 2; i++)
			guess = smoothers.get(level - 1)
			                 .mvMul(guess);
		Vector restrictedDefect =
			prolongationMatrices.get(level - 1)
			                    .tvMul(rhs.sub(systems.get(level)._1
				                                   .mvMul(guess)));
		final Vector correction = mgStep(level - 1,
		                                 new DenseVector(restrictedDefect.getLength()),
		                                 restrictedDefect);
		guess = guess.add(prolongationMatrices.get(level - 1)
		                                      .mvMul(correction));
		restrictedDefect =
			prolongationMatrices.get(level - 1)
			                    .tvMul(rhs.sub(systems.get(level)._1
				                                   .mvMul(guess)));
		
		for (int i = 0; i < 2; i++)
			guess = smoothers.get(level - 1)
			                 .mvMul(guess);
		//System.out.println("mgstep done" + level);
		return guess;
	}
	
	public Vector vCycle(final Vector initialIterate, final Vector rhs)
	{
		return mgStep(spaces.size() - 1, initialIterate, rhs);
	}
}
