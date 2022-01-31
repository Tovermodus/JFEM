package multigrid;

import basic.*;
import com.google.common.collect.TreeMultimap;
import io.vavr.Tuple2;
import linalg.*;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

public abstract class MGPreconditionerSpace<CSpace extends AcceptsMatrixBoundaryValues<CT, FT, ST, valueT, gradientT, hessianT
	> & Assembleable,
	CT extends Cell<CT, FT>,
	FT extends Face<CT, FT>,
	ST extends ShapeFunction<CT, FT, valueT, gradientT, hessianT> & Comparable<ST>, valueT, gradientT, hessianT>
	implements VectorMultiplyable
{
	public List<CSpace> spaces;
	public List<Smoother> smoothers;
	public List<VectorMultiplyable> prolongationOperators;
	public List<SparseMatrix> prolongationMatrices;
	public List<VectorMultiplyable> systems;
	public Vector finest_rhs;
	public VectorMultiplyable finest_system;
	public boolean verbose = false;
	public int cycles = 3;
	
	public MGPreconditionerSpace(final int refinements,
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
	
	VectorMultiplyable createProlongationMatrix(final CSpace coarse, final CSpace fine)
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
		final SparseMvMul prolong = new SparseMvMul(prolongationMatrix);
		
		return new VectorMultiplyable()
		{
			@Override
			public int getVectorSize()
			{
				return prolong.getVectorSize();
			}
			
			@Override
			public int getTVectorSize()
			{
				return prolong.getTVectorSize();
			}
			
			@Override
			public Vector mvMul(final Vector vector)
			{
				MutableVector v = new DenseVector(vector);
				applyZeroBoundaryConditions(coarse, v);
				v = prolongationMatrix.mvMul(v);
				return v;
			}
			
			@Override
			public Vector tvMul(final Vector vector)
			{
				final MutableVector v = prolong.tvMul(vector);
				applyZeroBoundaryConditions(coarse, v);
				return v;
			}
		};
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
	
	public abstract void applyCorrectBoundaryConditions(CSpace space, MutableVector vector);
	
	public abstract void applyZeroBoundaryConditions(CSpace space, MutableVector vector);
	
	public Vector mgStep(final int level, Vector guess, final Vector rhs)
	{
		if (level == 0)
		{
			if (verbose)
				System.out.println("       level 0 " + systems.get(0)
				                                              .getVectorSize());
			final Vector solution;
			if (systems.get(0) instanceof Matrix)
			{
				solution = new DenseMatrix((Matrix) systems.get(0)).solve(rhs);
			} else
			{
				solution = new IterativeSolver(true).solveGMRES(systems.get(0), rhs, 1e-9);
			}
			if (verbose)
				System.out.println("       solved" + solution.euclidianNorm());
			return solution;
		}
		if (verbose)
			System.out.println("   MG level " + level + " residual " + rhs.sub(systems.get(level)
			                                                                          .mvMul(guess))
			                                                              .euclidianNorm());
		guess = smoothers.get(level - 1)
		                 .smooth(systems.get(level), rhs, guess, verbose);
		Vector restrictedDefect =
			prolongationOperators.get(level - 1)
			                     .tvMul(rhs.sub(systems.get(level)
			                                           .mvMul(guess)));
		//restrictedDefect = applyCoarseBoundaryConditions(spaces.get(level - 1), restrictedDefect);
		if (verbose)
			System.out.println("   MG level " + level + " precorrection residual " + rhs.sub(systems.get(
				                                                                                        level)
			                                                                                        .mvMul(guess))
			                                                                            .euclidianNorm());
		if (verbose)
			System.out.println("   precorrection defect " + restrictedDefect.euclidianNorm());
		final Vector correction = mgStep(level - 1,
		                                 new DenseVector(restrictedDefect.getLength()),
		                                 restrictedDefect).mul(1);
		guess = guess.add(prolongationOperators.get(level - 1)
		                                       .mvMul(correction));
		if (verbose)
			restrictedDefect =
				prolongationOperators.get(level - 1)
				                     .tvMul(rhs.sub(systems.get(level)
				                                           .mvMul(guess)));
		if (verbose)
			System.out.println("    postcorrection defect " + restrictedDefect.euclidianNorm());
		if (verbose)
			System.out.println("   MG level " + level + " correction residual " + rhs.sub(systems.get(level)
			                                                                                     .mvMul(guess))
			                                                                         .euclidianNorm());
		guess = smoothers.get(level - 1)
		                 .smooth(systems.get(level), rhs, guess, verbose);
		if (verbose)
			System.out.println("   MG level " + level + " final residual " + rhs.sub(systems.get(level)
			                                                                                .mvMul(guess))
			                                                                    .euclidianNorm());
		return guess;
	}
	
	public Vector vCycle(final Vector initialIterate, final Vector rhs)
	{
		if (initialIterate.size() != rhs.size() || initialIterate.size() != getFinestSystem().getVectorSize())
			throw new IllegalArgumentException("wrong size");
		return mgStep(spaces.size() - 1, initialIterate, rhs);
	}
	
	public Vector fullVCycleCorrection(final Vector defect)
	{
		Vector d = new DenseVector(defect);
		final List<Vector> rhsides = new ArrayList<>();
		rhsides.add(d);
		for (int i = prolongationOperators.size() - 1; i >= 0; i--)
		{
			d = prolongationOperators.get(i)
			                         .tvMul(d);
			rhsides.add(d);
		}
		MutableVector iterate = new DenseVector(systems.get(0)
		                                               .getVectorSize());
		for (int i = 0; i < systems.size(); i++)
		{
			//applyZeroBoundaryConditions(spaces.get(i), iterate);
			iterate = new DenseVector(mgStep(i, iterate, rhsides.get(spaces.size() - i - 1)));
			//applyZeroBoundaryConditions(spaces.get(i), iterate);
			if (i < prolongationOperators.size())
				iterate = prolongationMatrices.get(i)
				                              .mvMul(iterate);
		}
		return iterate;
	}
	
	public VectorMultiplyable getFinestSystem()
	{
		return finest_system;
	}
	
	public CSpace getFinestSpace()
	{
		return spaces.get(spaces.size() - 1);
	}
	
	@Override
	public int getVectorSize()
	{
		return finest_system.getVectorSize();
	}
	
	@Override
	public int getTVectorSize()
	{
		return finest_system.getTVectorSize();
	}
	
	@Override
	public Vector mvMul(final Vector vector)
	{
		final MutableVector initial = new DenseVector(vector);
		final Vector defect = vector.sub(finest_system.mvMul(initial));
		System.out.println();
		Vector iterate = fullVCycleCorrection(defect);
		
		System.out.println("MG after FullVCycle " + finest_system.mvMul(initial.add(iterate))
		                                                         .sub(vector)
		                                                         .euclidianNorm() + " FROM " + defect.euclidianNorm());
		for (int i = 0; i < cycles; i++)
			iterate = vCycle(iterate, defect);
		
		System.out.println("MG after Second VCycle " + finest_system.mvMul(initial.add(iterate))
		                                                            .sub(vector)
		                                                            .euclidianNorm() + " FROM " + defect.euclidianNorm());
		return initial.add(iterate);
	}
	
	@Override
	public Vector tvMul(final Vector vector)
	{
		throw new UnsupportedOperationException("not implemented yet");
	}
}
