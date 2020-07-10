package multigrid;

import basic.MatrixFESpace;

public class RichardsonSmoother extends Smoother
{
//	private double omega;
//	public RichardsonSmoother(MatrixFESpace g, String[] args)
//	{
//		super(g, args);
//		this.omega = Double.parseDouble(args[0]);
//	}
//	@Override
//	public DoubleTensor smooth(DoubleTensor iterate, DoubleTensor rightHandSide)
//	{
//		DoubleTensor residuum = rightHandSide.sub(fESpace.getSystemMatrix().mvmul(iterate));
//		System.out.println(residuum.vectorNorm());
//		iterate = iterate.add(residuum).mul(omega);
//		return iterate;
//
//	}
}
