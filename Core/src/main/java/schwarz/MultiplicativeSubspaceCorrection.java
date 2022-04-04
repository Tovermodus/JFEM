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
			iterate = upgradePatch(schwarz, globalRhs, iterate, i);
		}
		for (int i = schwarz.getPatchCount() - 1; i >= 0; i--)
		{
			iterate = upgradePatch(schwarz, globalRhs, iterate, i);
		}
		return iterate;
	}
	
	private Vector upgradePatch(final AbstractSchwarz<?, ?, OT> schwarz,
	                            @NotNull final Vector globalRhs,
	                            Vector iterate,
	                            final int i)
	{
//		System.out.println();
//		System.out.println();
//		System.out.println();
		final Vector globalResidual = globalRhs.sub(schwarz.getGlobalOperator()
		                                                   .mvMul(iterate));
//		System.out.println(globalResidual.euclidianNorm());
		//System.out.println("SUB " + globalResidual.euclidianNorm());
		final Vector localRes = schwarz.getLocalVector(i,
		                                               globalResidual.mul(1));
//		System.out.println("res");
//		System.out.println(localRes.getLength());
//		System.out.println(globalResidual.getLength());
//		System.out.println(localRes.sub(globalResidual)
//		                           .euclidianNorm());
		final Vector localSol = schwarz.solveLocalSystem(i,
		                                                 localRes);
		final Vector globalSolComponent
			= schwarz.getGlobalVector(i,
			                          localSol);
//		System.out.println("sol");
//		System.out.println(localSol.sub(globalSolComponent)
//		                           .euclidianNorm());
//		System.out.println("acc");
//		System.out.println(schwarz.getLocalOperator(i)
//		                          .mvMul(localSol)
//		                          .sub(localRes)
//		                          .euclidianNorm());
//		System.out.println(schwarz.getGlobalOperator()
//		                          .mvMul(globalSolComponent)
//		                          .sub(globalResidual)
//		                          .euclidianNorm());
//		System.out.println("add");
//		System.out.println(globalRhs.sub(schwarz.getGlobalOperator()
//		                                        .mvMul(iterate.add(globalSolComponent)))
//		                            .euclidianNorm());
		if (i % 1 == 1 && localRes.euclidianNorm() > 1e-2)
		{
			final MixedTPFESpaceFunction<QkQkFunction> fun =
				new MixedTPFESpaceFunction<>(space.getShapeFunctionMap(),
				                             iterate);
			PlotWindow.addPlotShow(new MixedPlot2D(fun,
			                                       space.generatePlotPoints(32),
			                                       32));
		}
		iterate = iterate.add(globalSolComponent);
		return iterate;
	}
}
