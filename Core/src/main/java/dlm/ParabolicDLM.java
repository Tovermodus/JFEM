package dlm;

import basic.ScalarShapeFunction;
import basic.VectorShapeFunction;
import distorted.DistortedVectorCellIntegral;
import distorted.DistortedVectorSpace;
import linalg.BlockSparseMatrix;
import mixed.*;
import tensorproduct.CartesianGridSpace;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

import java.util.List;

public class ParabolicDLM<ST extends MixedShapeFunction<TPCell, TPFace, PF, VF>,
	PF extends ScalarShapeFunction<TPCell, TPFace>, VF extends VectorShapeFunction<TPCell, TPFace>>
{
	final CartesianGridSpace<ST, MixedValue, MixedGradient, MixedHessian> backgroundSpace;
	final List<DistortedVectorSpace> particleSpaces;
	BlockSparseMatrix constantSystemMatrix;
	BlockSparseMatrix currentSystemMatrix;
	List<MixedTPCellIntegral<PF, VF, ST>> constantBackgroundIntegrals;
	List<MixedTPCellIntegral<PF, VF, ST>> changingBackgroundIntegrals;
	List<List<DistortedVectorCellIntegral>> constantParticleIntegrals;
	List<List<DistortedVectorCellIntegral>> changingParticleIntegrals;
	
	public ParabolicDLM(final double dt,
	                    final int timeSteps,
	                    final CartesianGridSpace<ST, MixedValue, MixedGradient, MixedHessian> backgroundSpace,
	                    final List<DistortedVectorSpace> particleSpaces)
	{
		this.backgroundSpace = backgroundSpace;
		this.particleSpaces = particleSpaces;
	}
}
