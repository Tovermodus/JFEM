package basic;

import com.google.common.collect.Lists;
import io.vavr.Function3;

import java.util.List;
import java.util.function.BiFunction;
import java.util.function.DoubleSupplier;

public interface FESpaceTools<CT extends Cell<CT,FT>, FT extends  Face<CT,FT>,
	ST extends ShapeFunction<CT,FT,ST,?,?,?>>
{
	default void loopMatrixViaCell(Function3<CT,ST,ST, Double> integralEvaluation,
	                               MatrixFESpace<CT,FT,ST,?,?,?> space){
		List<List<CT>> smallerList = Lists.partition(space.getCells(),space.getCells().size()/12+1);
		smallerList.stream().parallel().forEach(smallList->
		{
			for (CT K : smallList)
			{
				for (ST v : space.getShapeFunctionsWithSupportOnCell(K))
				{
					for (ST u : space.getShapeFunctionsWithSupportOnCell(K))
					{
						double integral = integralEvaluation.apply(K,u,v);
						if(integral != 0)
							space.getSystemMatrix().add( integral,v.getGlobalIndex(),
								u.getGlobalIndex());
					}
				}
			}
		});
	}
	default void loopRhsViaCell(BiFunction<CT,ST,Double> integralEvaluation, MatrixFESpace<CT,FT,ST,?,?,?> space){
		List<List<CT>> smallerList = Lists.partition(space.getCells(),space.getCells().size()/12+1);
		smallerList.stream().parallel().forEach(smallList->
		{
			for (CT K : smallList)
			{
				for (ST v : space.getShapeFunctionsWithSupportOnCell(K))
				{
					double integral = integralEvaluation.apply(K,v);
					if(integral != 0)
						space.getRhs().add( integral,v.getGlobalIndex());
				}
			}
		});
	}
	default void loopMatrixViaFace(Function3<FT,ST,ST,Double> integralEvaluation,
	                               MatrixFESpace<CT,FT,ST,?,?,?> space){
		List<List<FT>> smallerList = Lists.partition(space.getFaces(),space.getFaces().size()/12+1);
		smallerList.stream().parallel().forEach(smallList->
		{
			for (FT F : smallList)
			{
				for (ST v : space.getShapeFunctionsWithSupportOnFace(F))
				{
					for (ST u : space.getShapeFunctionsWithSupportOnFace(F))
					{
						double integral = integralEvaluation.apply(F,u,v);
						if(integral != 0)
							space.getSystemMatrix().add( integral,v.getGlobalIndex(),
								u.getGlobalIndex());
					}
				}
			}
		});
	}
	default void loopRhsViaFace(BiFunction<FT,ST,Double> integralEvaluation, MatrixFESpace<CT,FT,ST,?,?,?> space){
		List<List<FT>> smallerList = Lists.partition(space.getFaces(),space.getFaces().size()/12+1);
		smallerList.stream().parallel().forEach(smallList->
		{
			for (FT F : smallList)
			{
				for (ST v : space.getShapeFunctionsWithSupportOnFace(F))
				{
					double integral = integralEvaluation.apply(F,v);
					if(integral != 0)
						space.getRhs().add( integral,v.getGlobalIndex());
				}
			}
		});
	}
}