package schwarz;

import basic.Cell;
import basic.Face;
import linalg.Vector;
import linalg.VectorMultiplyable;

import java.util.Collection;

public abstract class AbstractSchwarz<CT extends Cell<CT, FT>, FT extends Face<CT, FT>, OT extends VectorMultiplyable>
	implements VectorMultiplyable
{
	final private SubspaceCorrection<OT> subspaceCorrection;
	final private SystemSolver<OT> solver;
	
	protected AbstractSchwarz(final SubspaceCorrection<OT> subspaceCorrection, final SystemSolver<OT> solver)
	{
		this.subspaceCorrection = subspaceCorrection;
		this.solver = solver;
	}
	
	public abstract Collection<CT> getCellPatch(int patch);
	
	public abstract OT getGlobalOperator();
	
	public abstract OT getRestrictionOperator(int patch);
	
	public abstract OT getLocalOperator(int patch);
	
	public Vector getLocalVector(final int patch, final Vector globalVector)
	{
		return getRestrictionOperator(patch).mvMul(globalVector);
	}
	
	public Vector getGlobalVector(final int patch, final Vector globalVector)
	{
		return getRestrictionOperator(patch).tvMul(globalVector);
	}
	
	public abstract int getPatchCount();
	
	public SubspaceCorrection<OT> getSubspaceCorrection()
	{
		return subspaceCorrection;
	}
	
	public Vector solveLocalSystem(final int patch, final Vector localVector)
	{
		return solver.solve(getLocalOperator(patch), localVector);
	}
	
	@Override
	public int getVectorSize()
	{
		return getGlobalOperator().getVectorSize();
	}
	
	@Override
	public int getTVectorSize()
	{
		return getGlobalOperator().getTVectorSize();
	}
	
	@Override
	public Vector mvMul(final Vector vector)
	{
		return getSubspaceCorrection().solve(this, vector);
	}
	
	@Override
	public Vector tvMul(final Vector vector)
	{
		throw new IllegalArgumentException("Not implemented");
	}
}
