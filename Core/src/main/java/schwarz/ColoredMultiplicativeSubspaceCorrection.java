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
	
	public ColoredMultiplicativeSubspaceCorrection(final TaylorHoodSpace space)
	{
		this.space = space;
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
					                     if (j == 0)
						                     return globalResidual.mul(0);
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
				iterate = iterate.add(v);
		}
		return iterate;
	}
	
	@Override
	public Vector apply(final AbstractSchwarz<?, ?, OT> schwarz,
	                    @NotNull final Vector globalIterate,
	                    @NotNull final Vector globalRhs)
	{
		System.out.println("APPLY");
		Vector iterate = new DenseVector(globalIterate);
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
				iterate = iterate.add(v);
//				final MixedTPFESpaceFunction<QkQkFunction> fun =
//					new MixedTPFESpaceFunction<>(space.getShapeFunctionMap(), iterate);
//				PlotWindow.addPlotShow(new MixedPlot2D(fun,
//				                                       space,
//				                                       (int) Math.sqrt(space.getShapeFunctions()
//				                                                            .size()) / 2));
			}
		}
		return iterate;
	}
}
