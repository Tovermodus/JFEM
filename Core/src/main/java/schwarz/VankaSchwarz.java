package schwarz;

import it.unimi.dsi.fastutil.ints.IntOpenHashSet;
import it.unimi.dsi.fastutil.ints.IntSet;
import linalg.Matrix;
import linalg.SparseMvMul;
import mixed.TaylorHoodSpace;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

//pressure node oriented vanka
public class VankaSchwarz
	extends UpFrontSchwarz<TPCell, TPFace>
{
	final SubspaceCorrection<Matrix> subspaceCorrection;
	private final TaylorHoodSpace space;
	SparseMvMul globalSparsityPattern;
	
	public VankaSchwarz(final Matrix globalMatrix,
	                    final TaylorHoodSpace space,
	                    final SubspaceCorrection<Matrix> subspaceCorrection,
	                    final SystemSolver<Matrix> solver)
	{
		super(globalMatrix, subspaceCorrection, solver);
		this.space = space;
		this.subspaceCorrection = subspaceCorrection;
		if (space.getDimension() != 2)
			throw new IllegalArgumentException("Only in 2D");
		globalSparsityPattern = new SparseMvMul(globalMatrix);
		
		System.out.println("partitioned");
		build();
	}
	
	public IntSet getPatchDoFs(final int patch)
	{
		final int pressureDof = space.getVelocitySize() + patch;
		final int[] connectedDofs = globalSparsityPattern.getCol(pressureDof)._1;
		final IntSet s = new IntOpenHashSet();
		for (final int dof : connectedDofs)
			if (dof < space.getVelocitySize())
				s.add(dof);
		s.add(pressureDof);
		return s;
	}
	
	@Override
	public int getPatchCount()
	{
		return space.getShapeFunctions()
		            .size() - space.getVelocitySize();
	}
	
	@Override
	public SubspaceCorrection<Matrix> getSubspaceCorrection()
	{
		return subspaceCorrection;
	}
	
	@Override
	public RestrictionMatrix buildRestrictionMatrix(final int patch)
	{
		final IntSet patchDofs = getPatchDoFs(patch);
		final RestrictionMatrix s = new RestrictionMatrix(patchDofs.size(),
		                                                  space.getShapeFunctions()
		                                                       .size());
		int i = 0;
		for (final int dof : patchDofs)
			s.add(1, i++, dof);
		return s;
	}
}
