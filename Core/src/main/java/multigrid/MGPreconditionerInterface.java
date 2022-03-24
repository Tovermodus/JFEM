package multigrid;

import basic.*;
import com.google.common.collect.TreeMultimap;
import linalg.*;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.function.BiFunction;

public interface MGPreconditionerInterface<CSpace extends AcceptsMatrixBoundaryValues<CT, FT, ST, valueT, gradientT, hessianT
	> & Assembleable,
	CT extends Cell<CT, FT>,
	FT extends Face<CT, FT>,
	ST extends ShapeFunction<CT, FT, valueT, gradientT, hessianT> & Comparable<ST>, valueT, gradientT, hessianT>
	extends VectorMultiplyable
{
	CSpace getSpace(int level);
	
	Smoother getSmoother(int level);
	
	VectorMultiplyable getProlongationOperator(int level);
	
	SparseMatrix getProlongationMatrix(int level);
	
	VectorMultiplyable getSystem(int level);
	
	int maxLevel();
	
	boolean isVerbose();
	
	default int levelFromSize(final int size)
	{
		int level = maxLevel();
		while (getSystem(level)
			.getVectorSize() != size)
		{
			level--;
			if (level < 0)
				throw new IllegalArgumentException("defect has wrong shape");
		}
		return level;
	}
	
	default Vector restrictToSize(final int newSize, Vector finest)
	{
		int l = maxLevel() - 1;
		while (finest.getLength() != newSize)
		{
			if (l < 0)
				throw new IllegalArgumentException("newsize not valid");
			finest = new DenseVector(getProlongationOperator(l)
				                         .tvMul(finest));
			l--;
		}
		return finest;
	}
	
	default TreeMultimap<ST, ST> getRefinedFunctions(final CSpace coarse, final CSpace fine)
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
	
	default void presmoothcallback(final int level, final Vector guess, final Vector rhs)
	{
	
	}
	
	default void postmoothcallback(final int level, final Vector guess, final Vector rhs)
	{
	
	}
	
	default void correctioncallback(final int level, final Vector correction, final Vector rhs)
	{
	
	}
	
	default void precorrectioncallback(final int level, final Vector guess, final Vector rhs)
	{
	
	}
	
	default void postcorrectioncallback(final int level, final Vector guess, final Vector rhs)
	{
	
	}
	
	void applyZeroBoundaryConditions(CSpace space, MutableVector vector);
	
	default void applyBoundaryConditions(final CSpace space, final MutableVector vector,
	                                     final Vector boundaryConditionRhs)
	{
		applyZeroBoundaryConditions(space, vector);
		final DenseVector nonBoundaryRhs = new DenseVector(boundaryConditionRhs);
		applyZeroBoundaryConditions(space, nonBoundaryRhs);
		final Vector onlyBoundaryRhs = boundaryConditionRhs.sub(nonBoundaryRhs);
		vector.addInPlace(onlyBoundaryRhs);
	}
	
	default Vector mgStepV(final int level, Vector guess, final Vector rhs)
	{
		final String prefSpaces = "   .".repeat(maxLevel() - level + 1);
		if (level == 0)
		{
			if (isVerbose())
				System.out.println(prefSpaces + " level 0 " + getSystem(0)
					.getVectorSize());
			final Vector solution;
			if (getSystem(0) instanceof Matrix)
			{
				solution = new DenseMatrix((Matrix) getSystem(0)).solve(rhs);
			} else
			{
				solution = new IterativeSolver(true).solveGMRES(getSystem(0), rhs, 1e-9);
			}
			postmoothcallback(level, solution, rhs);
			if (isVerbose())
				System.out.println(prefSpaces + " solved " + getSystem(0).mvMul(solution)
				                                                         .sub(rhs)
				                                                         .euclidianNorm());
			return solution;
		}
		if (isVerbose())
			System.out.println(prefSpaces + " MG level " + level + " residual " + rhs.sub(getSystem(level)
				                                                                              .mvMul(guess))
			                                                                         .euclidianNorm());
		presmoothcallback(level, guess, rhs);
		guess = getSmoother(level)
			.smooth(getSystem(level), rhs, guess, isVerbose(), prefSpaces);
		precorrectioncallback(level, guess, rhs);
		Vector restrictedDefect =
			getProlongationOperator(level - 1)
				.tvMul(rhs.sub(getSystem(level)
					               .mvMul(guess)));
		//restrictedDefect = applyCoarseBoundaryConditions(spaces.get(level - 1), restrictedDefect);
		if (isVerbose())
			System.out.println(prefSpaces + " MG level " + level + " precorrection residual " + rhs.sub(
				                                                                                       getSystem(level)
					                                                                                       .mvMul(guess))
			                                                                                       .euclidianNorm());
		if (isVerbose())
			System.out.println(prefSpaces + " precorrection defect " + restrictedDefect.euclidianNorm());
		final Vector correction = mgStepV(level - 1,
		                                  new DenseVector(restrictedDefect.getLength()),
		                                  restrictedDefect);
		final Vector prolongedCorrection = getProlongationOperator(level - 1)
			.mvMul(correction);
		correctioncallback(level, prolongedCorrection, rhs);
		guess = guess.add(prolongedCorrection);
		postcorrectioncallback(level, guess, rhs);
		if (isVerbose())
			restrictedDefect =
				getProlongationOperator(level - 1)
					.tvMul(rhs.sub(getSystem(level)
						               .mvMul(guess)));
		if (isVerbose())
		{
			System.out.println(prefSpaces + " postcorrection defect " + restrictedDefect.euclidianNorm());
			System.out.println(prefSpaces + " MG level " + level + " correction residual " + rhs.sub(
				                                                                                    getSystem(level)
					                                                                                    .mvMul(guess))
			                                                                                    .euclidianNorm());
		}
		guess = getSmoother(level)
			.smooth(getSystem(level), rhs, guess, isVerbose(), prefSpaces);
		postmoothcallback(level, guess, rhs);
		if (isVerbose())
			System.out.println(prefSpaces + " MG level " + level + " final residual " + rhs.sub(getSystem(
				                                                                               level)
				                                                                                    .mvMul(guess))
			                                                                               .euclidianNorm());
		return guess;
	}
	
	default Vector mgStepW(final int level, Vector guess, final Vector rhs)
	{
		final String prefSpaces = "   .".repeat(maxLevel() - level + 1);
		if (level == 0)
		{
			if (isVerbose())
				System.out.println(prefSpaces + " level 0 " + getSystem(0)
					.getVectorSize());
			final Vector solution;
			if (getSystem(0) instanceof Matrix)
			{
				solution = new DenseMatrix((Matrix) getSystem(0)).solve(rhs);
			} else
			{
				solution = new IterativeSolver(true).solveGMRES(getSystem(0), rhs, 1e-9);
			}
			postmoothcallback(level, solution, rhs);
			if (isVerbose())
				System.out.println(prefSpaces + " solved" + solution.euclidianNorm());
			return solution;
		}
		if (isVerbose())
			System.out.println(prefSpaces + " MG level " + level + " residual " + rhs.sub(getSystem(level)
				                                                                              .mvMul(guess))
			                                                                         .euclidianNorm());
		presmoothcallback(level, guess, rhs);
		guess = getSmoother(level)
			.smooth(getSystem(level), rhs, guess, isVerbose(), prefSpaces);
		precorrectioncallback(level, guess, rhs);
		Vector restrictedDefect =
			getProlongationOperator(level - 1)
				.tvMul(rhs.sub(getSystem(level)
					               .mvMul(guess)));
		//restrictedDefect = applyCoarseBoundaryConditions(spaces.get(level - 1), restrictedDefect);
		if (isVerbose())
			System.out.println(prefSpaces + " MG level " + level + " precorrection residual " + rhs.sub(
				                                                                                       getSystem(level)
					                                                                                       .mvMul(guess))
			                                                                                       .euclidianNorm());
		if (isVerbose())
			System.out.println(prefSpaces + " precorrection defect " + restrictedDefect.euclidianNorm());
		Vector correction = mgStepV(level - 1,
		                            new DenseVector(restrictedDefect.getLength()),
		                            restrictedDefect).mul(1);
		Vector prolongedCorrection = getProlongationOperator(level - 1)
			.mvMul(correction);
		correctioncallback(level, prolongedCorrection, rhs);
		guess = guess.add(prolongedCorrection);
		postcorrectioncallback(level, guess, rhs);
		restrictedDefect =
			getProlongationOperator(level - 1)
				.tvMul(rhs.sub(getSystem(level)
					               .mvMul(guess)));
		//restrictedDefect = applyCoarseBoundaryConditions(spaces.get(level - 1), restrictedDefect);
		if (isVerbose())
			System.out.println(prefSpaces + " MG level " + level + " intercorrection residual " + rhs.sub(
				                                                                                         getSystem(level)
					                                                                                         .mvMul(guess))
			                                                                                         .euclidianNorm());
		if (isVerbose())
			System.out.println(prefSpaces + " intercorrection defect " + restrictedDefect.euclidianNorm());
		correction = mgStepV(level - 1,
		                     new DenseVector(restrictedDefect.getLength()),
		                     restrictedDefect).mul(1);
		prolongedCorrection = getProlongationOperator(level - 1)
			.mvMul(correction);
		correctioncallback(level, prolongedCorrection, rhs);
		guess = guess.add(prolongedCorrection);
		postcorrectioncallback(level, guess, rhs);
		if (isVerbose())
			restrictedDefect =
				getProlongationOperator(level - 1)
					.tvMul(rhs.sub(getSystem(level)
						               .mvMul(guess)));
		if (isVerbose())
		{
			System.out.println(prefSpaces + " postcorrection defect " + restrictedDefect.euclidianNorm());
			System.out.println(prefSpaces + " MG level " + level + " correction residual " + rhs.sub(
				                                                                                    getSystem(level)
					                                                                                    .mvMul(guess))
			                                                                                    .euclidianNorm());
		}
		guess = getSmoother(level)
			.smooth(getSystem(level), rhs, guess, isVerbose(), prefSpaces);
		postmoothcallback(level, guess, rhs);
		if (isVerbose())
			System.out.println(prefSpaces + " MG level " + level + " final residual " + rhs.sub(getSystem(
				                                                                               level)
				                                                                                    .mvMul(guess))
			                                                                               .euclidianNorm());
		return guess;
	}
	
	default Vector vCycle(final Vector initialIterate, final Vector rhs)
	{
		if (initialIterate.size() != rhs.size() || initialIterate.size() != getFinestSystem().getVectorSize())
			throw new IllegalArgumentException("wrong size");
		return mgStepV(maxLevel(), initialIterate, rhs);
	}
	
	default Vector wCycle(final Vector initialIterate, final Vector rhs)
	{
		if (initialIterate.size() != rhs.size() || initialIterate.size() != getFinestSystem().getVectorSize())
			throw new IllegalArgumentException("wrong size");
		return mgStepW(maxLevel(), initialIterate, rhs);
	}
	
	default Vector fullVCycleCorrection(final Vector defect)
	{
		Vector d = new DenseVector(defect);
		final List<Vector> rhsides = new ArrayList<>();
		rhsides.add(d);
		for (int i = maxLevel() - 1; i >= 0; i--)
		{
			d = getProlongationOperator(i)
				.tvMul(d);
			rhsides.add(d);
		}
		MutableVector iterate = new DenseVector(getSystem(0)
			                                        .getVectorSize());
		for (int i = 0; i < maxLevel(); i++)
		{
			iterate = new DenseVector(mgStepV(i, iterate, rhsides.get(maxLevel() - i)));
			iterate = getProlongationMatrix(i).mvMul(iterate);
		}
		return iterate;
	}
	
	default VectorMultiplyable getFinestSystem()
	{
		return getSystem(maxLevel());
	}
	
	default CSpace getFinestSpace()
	{
		return getSpace(maxLevel());
	}
	
	@Override
	default int getVectorSize()
	{
		return getFinestSystem().getVectorSize();
	}
	
	@Override
	default int getTVectorSize()
	{
		return getFinestSystem().getTVectorSize();
	}
	
	int getCycles();
	
	@Override
	default Vector mvMul(final Vector vector)
	{
		final MutableVector initial = new DenseVector(vector.getLength());
		getFinestSpace().copyBoundaryValues(vector, initial);
		final Vector defect = vector.sub(getFinestSystem().mvMul(initial));
		Vector iterate = fullVCycleCorrection(defect);
		if (isVerbose())
			System.out.println("MG after FullVCycle " + getFinestSystem().mvMul(initial.add(iterate))
			                                                             .sub(vector)
			                                                             .euclidianNorm() + " FROM " + defect.euclidianNorm());
		for (int i = 0; i < getCycles(); i++)
			iterate = vCycle(iterate, defect);
		
		if (isVerbose())
			System.out.println("MG after Second VCycle " + getFinestSystem().mvMul(initial.add(iterate))
			                                                                .sub(vector)
			                                                                .euclidianNorm() + " FROM " + defect.euclidianNorm());
		return initial.add(iterate);
	}
	
	@Override
	default Vector tvMul(final Vector vector)
	{
		throw new UnsupportedOperationException("not implemented yet");
	}
	
	default MGPreconditionerInterface<CSpace, CT, FT, ST, valueT, gradientT, hessianT> AMGFromMatrix(final SparseMatrix mat,
	                                                                                                 final BiFunction<Integer, SparseMatrix, Smoother> levelToSmoother,
	                                                                                                 final boolean isVerbose,
	                                                                                                 final int cycles)
	{
		final var me = this;
		final int level = me.levelFromSize(mat.getVectorSize());
		if (mat.getRows() != me.getSystem(level)
		                       .getVectorSize() || mat.getCols() != me.getSystem(level)
		                                                              .getTVectorSize())
			throw new IllegalArgumentException("Matrix needs to be of same shape as finest system");
		final List<SparseMatrix> reversedListOfSystems = new ArrayList<>();
		reversedListOfSystems.add(mat);
		//PlotWindow.addPlot(new MatrixPlot(finest_system, "finsys"));
		for (int i = 1; i <= level; i++)
		{
			final SparseMatrix coarser =
				new SparseMatrix(getProlongationMatrix(level - i)
					                 .tmMul(reversedListOfSystems.get(i - 1))
					                 .mmMul(getProlongationMatrix(level - i)));
			reversedListOfSystems.add(coarser);
		}
		final List<SparseMatrix> AMGSystems = new ArrayList<>();
		for (int i = level; i >= 0; i--)
			AMGSystems.add(reversedListOfSystems.get(i));
		final List<Smoother> smoothers = new ArrayList<>();
		for (int l = 1; l <= level; l++)
			smoothers.add(levelToSmoother.apply(l, AMGSystems.get(l)));
		
		return new MGPreconditionerInterface<CSpace, CT, FT, ST, valueT, gradientT, hessianT>()
		{
			@Override
			public CSpace getSpace(final int level)
			{
				return me.getSpace(level);
			}
			
			@Override
			public Smoother getSmoother(final int level)
			{
				return smoothers.get(level - 1);
			}
			
			@Override
			public VectorMultiplyable getProlongationOperator(final int level)
			{
				return me.getProlongationOperator(level);
			}
			
			@Override
			public SparseMatrix getProlongationMatrix(final int level)
			{
				return me.getProlongationMatrix(level);
			}
			
			@Override
			public SparseMatrix getSystem(final int level)
			{
				return AMGSystems.get(level);
			}
			
			@Override
			public int maxLevel()
			{
				return level;
			}
			
			@Override
			public boolean isVerbose()
			{
				return isVerbose;
			}
			
			@Override
			public void applyZeroBoundaryConditions(final CSpace space, final MutableVector vector)
			{
			}
			
			@Override
			public int getCycles()
			{
				return cycles;
			}
		};
	}
	
	default MGPreconditionerInterface<CSpace, CT, FT, ST, valueT, gradientT, hessianT> AMGFromVectorMultipliable(
		final VectorMultiplyable mat, final SparseMatrix coarsest,
		final BiFunction<Integer, VectorMultiplyable, Smoother> levelToSmoother, final boolean isVerbose)
	{
		if (mat.getVectorSize() != getVectorSize())
			throw new IllegalArgumentException("Sizes dont match" + mat.getVectorSize() + " != " + getVectorSize());
		if (coarsest.getCols() != getSystem(0).getVectorSize())
			throw new IllegalArgumentException("Coarsest Sizes dont match" + coarsest.getCols() +
				                                   " != " + getSystem(0).getVectorSize());
		
		final var me = this;
		
		final List<VectorMultiplyable> reversedListOfSystems = new ArrayList<>();
		reversedListOfSystems.add(mat);
		//PlotWindow.addPlot(new MatrixPlot(finest_system, "finsys"));
		for (int i = 1; i < maxLevel(); i++)
		{
			final VectorMultiplyable coarser =
				VectorMultiplyable.transposeConcatenate(getProlongationMatrix(maxLevel() - i),
				                                        VectorMultiplyable.concatenate(
					                                        reversedListOfSystems.get(i - 1),
					                                        getProlongationMatrix(maxLevel() - i)));
//
//			final SparseMatrix coarser =
//				new SparseMatrix(getProlongationMatrix(maxLevel() - i)
//					                 .tmMul(reversedListOfSystems.get(i - 1))
//					                 .mmMul(getProlongationMatrix(maxLevel() - i)));
			reversedListOfSystems.add(coarser);
		}
		reversedListOfSystems.add(coarsest);
		final List<VectorMultiplyable> AMGSystems = new ArrayList<>();
		for (int i = maxLevel(); i >= 0; i--)
			AMGSystems.add(reversedListOfSystems.get(i));
		
		return new MGPreconditionerInterface<CSpace, CT, FT, ST, valueT, gradientT, hessianT>()
		{
			@Override
			public CSpace getSpace(final int level)
			{
				return me.getSpace(level);
			}
			
			@Override
			public Smoother getSmoother(final int level)
			{
				return levelToSmoother.apply(level, getSystem(level));
			}
			
			@Override
			public VectorMultiplyable getProlongationOperator(final int level)
			{
				return me.getProlongationOperator(level);
			}
			
			@Override
			public SparseMatrix getProlongationMatrix(final int level)
			{
				return me.getProlongationMatrix(level);
			}
			
			@Override
			public VectorMultiplyable getSystem(final int level)
			{
				return AMGSystems.get(level);
			}
			
			@Override
			public int maxLevel()
			{
				return me.maxLevel();
			}
			
			@Override
			public boolean isVerbose()
			{
				return isVerbose;
			}
			
			@Override
			public void applyZeroBoundaryConditions(final CSpace space, final MutableVector vector)
			{
			}
			
			@Override
			public int getCycles()
			{
				return 3;
			}
		};
	}
}
