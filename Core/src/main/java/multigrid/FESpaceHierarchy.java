package multigrid;

import basic.*;
import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Iterables;
import com.google.common.collect.Multimap;
import basic.FESpace;

import java.lang.reflect.Constructor;
import java.util.ArrayList;
import java.util.Map;

public class FESpaceHierarchy
{
//	private ArrayList<FESpace> fESpaces;
//	private ArrayList<Multimap<Cell,Cell>> cellMaps;
//	private ArrayList<Multimap<Face, Face>> faceMaps;
//	private ArrayList<DoubleTensor> prolongationMatrices;
//	public FESpaceHierarchy(FESpace fESpace)
//	{
//		fESpaces = new ArrayList<>();
//		cellMaps = new ArrayList<>();
//		faceMaps = new ArrayList<>();
//		prolongationMatrices = new ArrayList<>();
//		getfESpaces().add(fESpace);
//	}
//	public void addGloballyRefinedLevel()
//	{
//		FESpace lastFESpace = Iterables.getLast(getfESpaces());
//		Multimap<Cell,Cell> cellMap = ArrayListMultimap.create();
//		Multimap<Face,Face> faceMap = ArrayListMultimap.create();
//		FESpace refinedFESpace = lastFESpace.refineAll(cellMap,faceMap);
//		getfESpaces().add(refinedFESpace);
//		getCellMaps().add(cellMap);
//		getFaceMaps().add(faceMap);
//		DoubleTensor prolongationMatrix = new DoubleTensor(refinedFESpace.getShapeFunctions().size(),
//			lastFESpace.getShapeFunctions().size(),true);
//		getProlongationMatrices().add(prolongationMatrix);
//		for(ScalarShapeFunction shapeFunction: lastFESpace.getShapeFunctions())
//		{
//			ArrayList<ScalarShapeFunction> refinedFunctions = new ArrayList<>();
//			for(Cell cell:shapeFunction.getCells())
//			{
//				for(Cell refinedSupportCell: cellMap.get(cell))
//				{
//					refinedFunctions.addAll(refinedSupportCell.getShapeFunctions());
//				}
//			}
//			Map<Integer, Double> coeffs = shapeFunction.prolongate(refinedFunctions);
//			for(int i:coeffs.keySet())
//			{
//				prolongationMatrix.add(i, shapeFunction.getGlobalIndex(),coeffs.get(i));
//			}
//		}
//	}
//	public void evaluateCellIntegrals(ArrayList<CellIntegral> cellIntegrals, ArrayList<RightHandSideIntegral> rightHandSideIntegrals)
//	{
//		for(FESpace g: getfESpaces())
//			g.evaluateCellIntegrals(cellIntegrals,rightHandSideIntegrals);
//	}
//	public void evaluateFaceIntegrals(ArrayList<FaceIntegral> faceIntegrals, ArrayList<BoundaryFaceIntegral> boundaryFaceIntegrals)
//	{
//		for(FESpace g: getfESpaces())
//			g.evaluateFaceIntegrals(faceIntegrals,boundaryFaceIntegrals);
//	}
//
//	public DoubleTensor multiGridSolve(double tolerance, int maxiter,
//	                                                       Class<? extends Smoother>  smoother,String[] args)
//	{
//		ArrayList<Smoother> smoothers = new ArrayList<>();
//		Constructor<? extends Smoother> constr = null;
//		try
//		{
//			constr = smoother.getConstructor(FESpace.class, String[].class);
//		} catch (NoSuchMethodException e)
//		{
//			e.printStackTrace();
//		}
//		for (FESpace fESpace : getfESpaces())
//		{
//			try
//			{
//				assert constr != null;
//				smoothers.add(constr.newInstance(fESpace, args));
//				System.out.println("kjh");
//			} catch (Exception e)
//			{
//				e.printStackTrace();
//			}
//		}
//		FESpace lastFESpace = Iterables.getLast(getfESpaces());
//		System.out.println(lastFESpace.getSystemMatrix().getM()+" "+ lastFESpace.getSystemMatrix().getN());
//		DoubleTensor iterate = new DoubleTensor((int) lastFESpace.getRhs().size());
//		DoubleTensor res = new DoubleTensor(lastFESpace.getRhs());
//		DoubleTensor correct;
//		System.out.println(res.vectorNorm());
//		//res.print_formatted();
//		for (int i = 0; i < maxiter && res.vectorNorm() > tolerance; i++)
//		{
//			correct = multiGridV(getProlongationMatrices().size(), res, smoothers);
//			//correct.print_formatted("correct");
//			iterate = iterate.add(correct);
//			//iterate.print_formatted("mgiterate");
//			res = lastFESpace.getRhs().sub(lastFESpace.getSystemMatrix().mvmul(iterate));
//			//res.print_formatted("res");
//			System.out.println(i + " " + res.vectorNorm());
//		}
//		return iterate;
//	}
//	public DoubleTensor multiGridV(int level, DoubleTensor rightHandSide, ArrayList<Smoother> smoothers)
//	{
//		int preIters = 2;
//		int postIters = 2;
//		DoubleTensor iterate = new DoubleTensor(rightHandSide.size());
//		FESpace g = getfESpaces().get(level);
//		Smoother smoother = smoothers.get(level);
//		//rightHandSide.print_formatted("rhs");
//		if(level == 0)
//		{
//			//g.A.solve(rightHandSide).print_formatted("level0solve");
//			return g.getSystemMatrix().solve(rightHandSide);
//		}
//		for(int k = 0; k < preIters; k++)
//			iterate = smoother.smooth(iterate,rightHandSide);
//		DoubleTensor restrictedRightHandSide =
//			getProlongationMatrices().get(level - 1).tvmul(rightHandSide.sub(g.getSystemMatrix().mvmul(iterate)));
//		DoubleTensor correction = multiGridV(level - 1,
//			restrictedRightHandSide,smoothers);
//		iterate = iterate.add(getProlongationMatrices().get(level - 1).mvmul(correction));
//		//System.out.println("kajsdka");
//		//System.out.println(rightHandSide.sub(g.A.mvmul(iterate)).vectorNorm());
//		for(int k = 0; k < postIters; k++)
//			iterate = smoother.smooth(iterate,rightHandSide);
//		//System.out.println(rightHandSide.sub(g.A.mvmul(iterate)).vectorNorm());
//		//iterate.print_formatted("return");
//		return iterate;
//
//	}
//
//	public ArrayList<FESpace> getfESpaces()
//	{
//		return fESpaces;
//	}
//
//	public ArrayList<Multimap<Cell, Cell>> getCellMaps()
//	{
//		return cellMaps;
//	}
//
//	public ArrayList<Multimap<Face, Face>> getFaceMaps()
//	{
//		return faceMaps;
//	}
//
//	public ArrayList<DoubleTensor> getProlongationMatrices()
//	{
//		return prolongationMatrices;
//	}
}
