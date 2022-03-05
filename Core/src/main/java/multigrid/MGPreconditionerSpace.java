package multigrid;

import basic.*;
import com.google.common.collect.TreeMultimap;
import io.vavr.Tuple2;
import linalg.*;

import java.util.ArrayList;
import java.util.List;

public abstract class MGPreconditionerSpace<CSpace extends AcceptsMatrixBoundaryValues<CT, FT, ST, valueT, gradientT, hessianT
	> & Assembleable,
	CT extends Cell<CT, FT>,
	FT extends Face<CT, FT>,
	ST extends ShapeFunction<CT, FT, valueT, gradientT, hessianT> & Comparable<ST>, valueT, gradientT, hessianT>
	implements MGPreconditionerInterface<CSpace, CT, FT, ST, valueT, gradientT, hessianT>
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
			System.out.println("mg" + ((SparseMatrix) mat_rhs._1).sub(((SparseMatrix) mat_rhs._1).transpose())
			                                                     .absMaxElement());
			finest_rhs = mat_rhs._2;
			finest_system = mat_rhs._1;
			System.out.println("functions, system " + k++);
		}
		smoothers = createSmoothers();
		System.out.println("smoothers");
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
		//System.out.println(prolongationMatrix.mul(64));
		//PlotWindow.addPlot(new MatrixPlot(prolongationMatrix));
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
				v = prolongationMatrix.mvMul(v);
				applyZeroBoundaryConditions(fine, v);
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
	
	public abstract List<CSpace> createSpaces(int refinements);
	
	public abstract Tuple2<VectorMultiplyable, DenseVector> createSystem(CSpace space);
	
	public abstract List<Smoother> createSmoothers();
	
	@Override
	public CSpace getSpace(final int level)
	{
		return spaces.get(level);
	}
	
	@Override
	public Smoother getSmoother(final int level)
	{
		return smoothers.get(level - 1);
	}
	
	@Override
	public VectorMultiplyable getProlongationOperator(final int level)
	{
		return prolongationOperators.get(level);
	}
	
	@Override
	public SparseMatrix getProlongationMatrix(final int level)
	{
		return prolongationMatrices.get(level);
	}
	
	@Override
	public VectorMultiplyable getSystem(final int level)
	{
		return systems.get(level);
	}
	
	@Override
	public int maxLevel()
	{
		return spaces.size() - 1;
	}
	
	@Override
	public boolean isVerbose()
	{
		return verbose;
	}
	
	@Override
	public int getCycles()
	{
		return cycles;
	}
	
	Vector mvMulW(final Vector vector)
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
			iterate = wCycle(iterate, defect);
		
		if (isVerbose())
			System.out.println("MG after Second VCycle " + getFinestSystem().mvMul(initial.add(iterate))
			                                                                .sub(vector)
			                                                                .euclidianNorm() + " FROM " + defect.euclidianNorm());
		return initial.add(iterate);
	}
}
