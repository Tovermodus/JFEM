package basic;

import com.google.common.collect.Lists;
import linalg.Matrix;
import linalg.Vector;

import java.util.Arrays;
import java.util.List;

public interface MatrixFESpace<CT extends Cell<CT,FT,ST>, FT extends  Face<CT,FT,ST>,
	ST extends ShapeFunction<CT,FT,ST,valueT,gradientT,hessianT>,valueT,gradientT,hessianT,
	FST extends MatrixFESpace<CT,FT,ST,valueT,gradientT,hessianT,FST>> extends FESpace<CT,FT
	,ST,valueT,gradientT,hessianT,FST>
{
	void initializeSystemMatrix();
	void initializeRhs();
	
	Vector getRhs();
	
	Matrix getSystemMatrix();
	
	default void evaluateCellIntegrals(List<CellIntegral<CT,FT,ST>> cellIntegrals,
	                           List<RightHandSideIntegral<CT,FT,ST>> rightHandSideIntegrals)
		{
		List<List<CT>> smallerList = Lists.partition(getCells(),12);
		smallerList.stream().parallel().forEach(smallList->
		{
			for (CT K : smallList)
			{
				for (ST v : K.getShapeFunctions())
				{
					for (ST u : K.getShapeFunctions())
					{
						double integral = 0;
						for (CellIntegral<CT,FT,ST> cellIntegral :
							cellIntegrals)
						{
							integral += cellIntegral.evaluateCellIntegral(K, u, v);
						}
						if(integral != 0)
							getSystemMatrix().add( integral,v.getGlobalIndex(),
								u.getGlobalIndex());
					}
					double integral = 0;
					for (RightHandSideIntegral<CT,FT,ST> rightHandSideIntegral :
						rightHandSideIntegrals)
					{
						integral += rightHandSideIntegral.evaluateRightHandSideIntegral(K, v);
					}
					if(integral != 0)
						getRhs().add(integral, v.getGlobalIndex());

				}
			}

		});
	}
	default void evaluateFaceIntegrals(List<FaceIntegral<CT,FT,ST>> faceIntegrals,
	                                   List<BoundaryRightHandSideIntegral<CT,FT,ST>> boundaryRightHandSideIntegrals)
	{
		
		List<List<FT>> smallerList = Lists.partition(getFaces(),12);
//		for (ST s: getShapeFunctions())
//		{
//			if(s.getGlobalIndex() == 2)
//				s.setGlobalIndex(4);
//			else if(s.getGlobalIndex() == 3)
//				s.setGlobalIndex(5);
//			else if(s.getGlobalIndex() == 4)
//				s.setGlobalIndex(2);
//			else if(s.getGlobalIndex() == 5)
//				s.setGlobalIndex(3);
//			else if(s.getGlobalIndex() == 10)
//				s.setGlobalIndex(12);
//			else if(s.getGlobalIndex() == 11)
//				s.setGlobalIndex(13);
//			else if(s.getGlobalIndex() == 12)
//				s.setGlobalIndex(10);
//			else if(s.getGlobalIndex() == 13)
//				s.setGlobalIndex(11);
//		}
		smallerList.stream().parallel().forEach(smallList->
		{
			for (FT F : smallList)
			{
				for (ST u : F.getShapeFunctions())
				{
					
					for (ST v : F.getShapeFunctions())
					{
						double integral = 0;
						for (FaceIntegral<CT,FT,ST> faceIntegral :
							faceIntegrals)
						{
							integral += faceIntegral.evaluateFaceIntegral(F,u,v);
						}
						if(integral != 0)
							getSystemMatrix().add(integral,v.getGlobalIndex(),
								u.getGlobalIndex());
					}
					double integral = 0;
					for (BoundaryRightHandSideIntegral<CT,FT,ST> rightHandSideIntegral :
						boundaryRightHandSideIntegrals)
					{
						integral += rightHandSideIntegral.evaluateBoundaryRightHandSideIntegral(F, u);
					}
					if(integral != 0)
						getRhs().add(integral, u.getGlobalIndex());
				}
			}
		});
	}
}
