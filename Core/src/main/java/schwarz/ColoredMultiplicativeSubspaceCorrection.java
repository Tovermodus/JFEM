package schwarz;

import it.unimi.dsi.fastutil.ints.IntSet;
import linalg.DenseVector;
import linalg.Vector;
import linalg.VectorMultiplyable;
import mixed.TaylorHoodSpace;
import org.jetbrains.annotations.NotNull;

import java.util.List;
import java.util.stream.Collectors;

public class ColoredMultiplicativeSubspaceCorrection<OT extends VectorMultiplyable>
	implements SubspaceCorrection<OT>
{
	List<IntSet> colorSets;
	TaylorHoodSpace space;
	final double omega;
	
	public ColoredMultiplicativeSubspaceCorrection(final TaylorHoodSpace space, final double omega)
	{
		this.space = space;
		this.omega = omega;
	}
	
	public void setColors(final List<IntSet> colorSets)
	{
		
		this.colorSets = colorSets;
	}
	
	@Override
	public Vector solve(final AbstractSchwarz<?, ?, OT> schwarz, final Vector globalRhs)
	{
		Vector iterate = schwarz.getGlobalVector(0, schwarz.solveLocalSystem(0,
		                                                                     schwarz.getLocalVector(0,
		                                                                                            globalRhs)));
		for (int i = 0; i < colorSets.size(); i++)
		{
			final Vector globalResidual = globalRhs.sub(schwarz.getGlobalOperator()
			                                                   .mvMul(iterate));
			final IntSet coloredPatches = colorSets.get(i);
			final List<Vector> globalSolComponents
				= coloredPatches.stream()
				                .parallel()
				                .map(j ->
				                     {
//					                     if (j == 0)
//						                     return globalResidual.mul(0);
					                     final Vector localRes = schwarz.getLocalVector(j,
					                                                                    globalResidual);
					                     final Vector localSol = schwarz.solveLocalSystem(j,
					                                                                      localRes);
					                     final Vector globalSolComponent
						                     = schwarz.getGlobalVector(j,
						                                               localSol);
					                     return globalSolComponent;
				                     })
				                .collect(Collectors.toList());
			for (final Vector v : globalSolComponents)
				iterate = iterate.add(v.mul(omega));
		}
		return iterate;
	}
	
	@Override
	public Vector apply(final AbstractSchwarz<?, ?, OT> schwarz,
	                    @NotNull final Vector globalIterate,
	                    @NotNull final Vector globalRhs)
	{
		Vector iterate = new DenseVector(globalIterate);
		for (int i = 0; i < colorSets.size(); i++)
		{
			iterate = colorSetUpgrade((AbstractSchwarz<?, ?, OT>) schwarz, globalRhs, iterate, i);
		}
		for (int i = colorSets.size() - 1; i >= 0; i--)
		{
			iterate = colorSetUpgrade((AbstractSchwarz<?, ?, OT>) schwarz, globalRhs, iterate, i);
		}
		return iterate;
	}
	
	private Vector colorSetUpgrade(final AbstractSchwarz<?, ?, OT> schwarz,
	                               @NotNull final Vector globalRhs,
	                               Vector iterate,
	                               final int i)
	{
		final Vector globalResidual = globalRhs.sub(schwarz.getGlobalOperator()
		                                                   .mvMul(iterate));
		final IntSet coloredPatches = colorSets.get(i);
		final List<Vector> globalSolComponents
			= coloredPatches.stream()
			                .parallel()
			                .map(j ->
			                     {
				                     final Vector localRes
					                     = schwarz.getLocalVector(
					                     j,
					                     globalResidual);
				                     final Vector localSol
					                     = schwarz.solveLocalSystem(
					                     j,
					                     localRes);
				                     final Vector
					                     globalSolComponent
					                     = schwarz.getGlobalVector(
					                     j,
					                     localSol);
				                     return globalSolComponent;
			                     })
			                .collect(Collectors.toList());
		for (final Vector v : globalSolComponents)
		{
			iterate = iterate.add(v.mul(omega));
//				final MixedTPFESpaceFunction<QkQkFunction> fun =
//					new MixedTPFESpaceFunction<>(space.getShapeFunctionMap(), iterate);
//				PlotWindow.addPlotShow(new MixedPlot2D(fun,
//				                                       space,
//				                                       (int) Math.sqrt(space.getShapeFunctions()
//				                                                            .size()) / 2));
		}
		return iterate;
	}
}
