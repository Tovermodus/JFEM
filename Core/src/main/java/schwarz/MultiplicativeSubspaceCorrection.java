package schwarz;

import basic.PlotWindow;
import linalg.DenseVector;
import linalg.Vector;
import linalg.VectorMultiplyable;
import mixed.MixedPlot2D;
import mixed.MixedTPFESpaceFunction;
import mixed.QkQkFunction;
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
			if (i % 1 == 1 && localRes.euclidianNorm() > 1e-2)
			{
				final MixedTPFESpaceFunction<QkQkFunction> fun =
					new MixedTPFESpaceFunction<>(space.getShapeFunctionMap(), iterate);
				PlotWindow.addPlotShow(new MixedPlot2D(fun,
				                                       space.generatePlotPoints(20),
				                                       20));
			}
			iterate = iterate.add(globalSolComponent);
		}
		return iterate;
	}
}
