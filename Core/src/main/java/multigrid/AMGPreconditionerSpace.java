package multigrid;

import basic.*;
import com.google.common.collect.TreeMultimap;
import io.vavr.Tuple2;
import linalg.*;
import org.jetbrains.annotations.NotNull;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

public abstract class AMGPreconditionerSpace<CSpace extends AcceptsMatrixBoundaryValues<CT, FT, ST, valueT, gradientT, hessianT
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
	public List<SparseMatrix> systems;
	public Vector finest_rhs;
	public SparseMatrix finest_system;
	public boolean verbose = false;
	public int cycles = 3;
	
	public AMGPreconditionerSpace(final int refinements,
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
		final var mat_rhs = createSystem(spaces.get(spaces.size() - 1));
		finest_rhs = mat_rhs._2;
		finest_system = mat_rhs._1;
		System.out.println("finest sys");
		systems = createSystems();
		smoothers = createSmoothers();
		System.out.println("smoothers");
	}
	
	private List<SparseMatrix> createSystems()
	{
		final List<SparseMatrix> reversed = new ArrayList<>();
		reversed.add(finest_system);
		//PlotWindow.addPlot(new MatrixPlot(finest_system, "finsys"));
		for (int i = 1; i < spaces.size(); i++)
		{
			System.out.println("system " + i);
			final int systemIndex = spaces.size() - i - 1;
			final SparseMatrix coarser =
				new SparseMatrix(prolongationMatrices.get(systemIndex)
				                                     .tmMul(reversed.get(i - 1))
				                                     .mmMul(prolongationMatrices.get(systemIndex)));
			reversed.add(coarser);
			//PlotWindow.addPlot(new MatrixPlot(coarser, "sys " + i + " tem"));
			PlotWindow.addPlot(new MatrixPlot(prolongationMatrices.get(i - 1),
			                                  "prol " + i + " tem " + prolongationMatrices.get(i - 1)
			                                                                              .getCols()));
		}
		final List<SparseMatrix> ret = new ArrayList<>();
		for (int i = spaces.size() - 1; i >= 0; i--)
			ret.add(reversed.get(i));
		return ret;
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
		final SparseMatrix prolongationMatrix = buildProlongationMatrix(coarse, fine);
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
	
	@NotNull
	protected SparseMatrix buildProlongationMatrix(final CSpace coarse, final CSpace fine)
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
		return prolongationMatrix;
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
						                                          .get(coarseCell);
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
	
	public abstract Tuple2<SparseMatrix, DenseVector> createSystem(CSpace space);
	
	public abstract List<Smoother> createSmoothers();
	
	public abstract void applyZeroBoundaryConditions(CSpace space, MutableVector vector);
	
	public Vector mgStep(final int level, Vector guess, final Vector rhs)
	{
		final String prefSpaces = "   .".repeat(spaces.size() - level);
		if (level == 0)
		{
			if (verbose)
				System.out.println(prefSpaces + "level 0 " + systems.get(0)
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
				System.out.println(prefSpaces + "solved" + solution.euclidianNorm());
			return solution;
		}
		if (verbose)
			System.out.println(prefSpaces + "MG level " + level + " residual " + rhs.sub(systems.get(level)
			                                                                                    .mvMul(guess))
			                                                                        .euclidianNorm());
		guess = smoothers.get(level - 1)
		                 .smooth(systems.get(level), rhs, guess, verbose, prefSpaces);
		Vector restrictedDefect =
			prolongationMatrices.get(level - 1)
			                    .tvMul(rhs.sub(systems.get(level)
			                                          .mvMul(guess)));
		//restrictedDefect = applyCoarseBoundaryConditions(spaces.get(level - 1), restrictedDefect);
		if (verbose)
			System.out.println(prefSpaces + "MG level " + level + " precorrection residual " + rhs.sub(
				                                                                                      systems.get(
					                                                                                             level)
				                                                                                             .mvMul(guess))
			                                                                                      .euclidianNorm());
		if (verbose)
			System.out.println(prefSpaces + "precorrection defect " + restrictedDefect.euclidianNorm());
		final Vector correction = mgStep(level - 1,
		                                 new DenseVector(restrictedDefect.getLength()),
		                                 restrictedDefect).mul(1);
		guess = guess.add(prolongationMatrices.get(level - 1)
		                                      .mvMul(correction));
		if (verbose)
			restrictedDefect =
				prolongationMatrices.get(level - 1)
				                    .tvMul(rhs.sub(systems.get(level)
				                                          .mvMul(guess)));
		if (verbose)
			System.out.println(prefSpaces + "postcorrection defect " + restrictedDefect.euclidianNorm());
		if (verbose)
			System.out.println(prefSpaces + "MG level " + level + " correction residual " + rhs.sub(systems.get(
				                                                                                               level)
			                                                                                               .mvMul(guess))
			                                                                                   .euclidianNorm());
		guess = smoothers.get(level - 1)
		                 .smooth(systems.get(level), rhs, guess, verbose, prefSpaces);
		if (verbose)
			System.out.println(prefSpaces + "MG level " + level + " final residual " + rhs.sub(systems.get(
				                                                                                          level)
			                                                                                          .mvMul(guess))
			                                                                              .euclidianNorm());
		return guess;
	}
	
	public int levelFromSize(final int size)
	{
		int level = systems.size() - 1;
		while (systems.get(level)
		              .getCols() != size)
		{
			level--;
			if (level < 0)
				throw new IllegalArgumentException("defect has wrong shape");
		}
		return level;
	}
	
	public Vector vCycle(final Vector initialIterate, final Vector rhs)
	{
		if (initialIterate.size() != rhs.size())
			throw new IllegalArgumentException("wrong size");
		return mgStep(levelFromSize(initialIterate.getLength()), initialIterate, rhs);
	}
	
	public Vector fullVCycleCorrection(final Vector defect)
	{
		final int level = levelFromSize(defect.getLength());
		Vector d = new DenseVector(defect);
		final List<Vector> rhsides = new ArrayList<>();
		rhsides.add(d);
		for (int i = level - 1; i >= 0; i--)
		{
			d = prolongationOperators.get(i)
			                         .tvMul(d);
			rhsides.add(d);
		}
		MutableVector iterate = new DenseVector(systems.get(0)
		                                               .getVectorSize());
		for (int i = 0; i < level; i++)
		{
			iterate = new DenseVector(mgStep(i, iterate, rhsides.get(rhsides.size() - i - 1)));
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
		final int level = levelFromSize(vector.getLength());
		final MutableVector vectorWith0Boundary = new DenseVector(vector);
		applyZeroBoundaryConditions(spaces.get(level), vectorWith0Boundary);
		final Vector initial = vector.sub(vectorWith0Boundary);
		final Vector defect = vector.sub(systems.get(level)
		                                        .mvMul(initial));
		Vector iterate = fullVCycleCorrection(defect);
		if (verbose)
			System.out.println("MG after FullVCycle " + systems.get(level)
			                                                   .mvMul(initial.add(iterate))
			                                                   .sub(vector)
			                                                   .euclidianNorm() + " FROM " + defect.euclidianNorm());
		for (int i = 0; i < cycles; i++)
			iterate = vCycle(iterate, defect);
		
		if (verbose)
			System.out.println("MG after Second VCycle " + systems.get(level)
			                                                      .mvMul(initial.add(iterate))
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
