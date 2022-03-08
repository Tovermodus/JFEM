package schwarz;

import linalg.DenseVector;
import linalg.Vector;
import linalg.VectorMultiplyable;
import mixed.TaylorHoodSpace;
import org.jetbrains.annotations.NotNull;

public class MultiplicativeSubspaceCorrection<OT extends VectorMultiplyable>
	implements SubspaceCorrection<OT>
{
	TaylorHoodSpace space;
	
	public MultiplicativeSubspaceCorrection()
	{
	
	}
	
	public MultiplicativeSubspaceCorrection(final TaylorHoodSpace space)
	{
		this.space = space;
	}
	
	@Override
	public Vector solve(final AbstractSchwarz<?, ?, OT> schwarz, final Vector globalRhs)
	{
		Vector iterate = schwarz.getGlobalVector(0, schwarz.solveLocalSystem(0,
		                                                                     schwarz.getLocalVector(0,
		                                                                                            globalRhs)));
		for (int i = 1; i < schwarz.getPatchCount(); i++)
		{
			
			final Vector globalResidual = globalRhs.sub(schwarz.getGlobalOperator()
			                                                   .mvMul(iterate));
			final Vector localRes = schwarz.getLocalVector(i,
			                                               globalResidual);
			final Vector localSol = schwarz.solveLocalSystem(i,
			                                                 localRes);
			final Vector globalSolComponent
				= schwarz.getGlobalVector(i,
				                          localSol);
			iterate = iterate.add(globalSolComponent);
		}
		return iterate;
	}
	
	@Override
	public Vector apply(final AbstractSchwarz<?, ?, OT> schwarz,
	                    @NotNull final Vector globalIterate,
	                    @NotNull final Vector globalRhs)
	{
		Vector iterate = new DenseVector(globalIterate);
		final var points = space.generatePlotPoints((int) Math.sqrt(space.getShapeFunctions()
		                                                                 .size()) / 2);
		for (int i = 0; i < schwarz.getPatchCount(); i++)
		{
			final Vector globalResidual = globalRhs.sub(schwarz.getGlobalOperator()
			                                                   .mvMul(iterate));
			//System.out.println("SUB " + globalResidual.euclidianNorm());
			final Vector localRes = schwarz.getLocalVector(i,
			                                               globalResidual);
			final Vector localSol = schwarz.solveLocalSystem(i,
			                                                 localRes);
			final Vector globalSolComponent
				= schwarz.getGlobalVector(i,
				                          localSol);
//			if (i % 4 == 0)
//			{
//			final MixedTPFESpaceFunction<QkQkFunction> fun =
//				new MixedTPFESpaceFunction<>(space.getShapeFunctionMap(), iterate);
//			PlotWindow.addPlotShow(new MixedPlot2D(fun,
//			                                       points,
//			                                       (int) Math.sqrt(space.getShapeFunctions()
//			                                                            .size()) / 2));
//			}
			iterate = iterate.add(globalSolComponent);
		}
		return iterate;
	}
}
