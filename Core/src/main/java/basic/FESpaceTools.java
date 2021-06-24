package basic;

import com.google.common.collect.Lists;
import io.vavr.Function3;

import java.util.List;
import java.util.function.BiFunction;
import java.util.stream.Stream;

public interface FESpaceTools<CT extends Cell<CT,FT,ET>, FT extends  Face<CT,FT,ET>, ET extends Edge<CT,FT,ET>,
	ST extends ShapeFunction<CT,FT,ET,ST,?,?,?>>
{
	default int increaseCellCounter()
	{
		return 0;
	}
	default int increaseFaceCounter()
	{
		return 0;
	}
	default void loopMatrixViaCell(Function3<CT,ST,ST, Double> integralEvaluation,
	                               MatrixFESpace<CT,FT,ET,ST,?,?,?> space){
		List<List<CT>> smallerList = Lists.partition(space.getCells(),
			space.getCells().size()/PerformanceArguments.getInstance().threadNumber+1);
		Stream<List<CT>> stream = smallerList.stream();
		if(PerformanceArguments.getInstance().parallelizeThreads)
			stream = stream.parallel();
		stream.forEach(smallList->
		{
			for (CT K : smallList)
			{
				System.out.println((int)(100.0*increaseCellCounter()/space.getCells().size())+"%");
				for (ST v : space.getShapeFunctionsWithSupportOnCell(K))
				{
					if(space.getFixedNodeIndices().contains(v.getGlobalIndex()))
						continue;
					for (ST u : space.getShapeFunctionsWithSupportOnCell(K))
					{
						if(space.getFixedNodeIndices().contains(u.getGlobalIndex()))
							continue;
						if(space.getFixedNodeIndices().contains(u.getGlobalIndex()))
							continue;
						double integral = integralEvaluation.apply(K,u,v);
						if(integral != 0)
							space.getSystemMatrix().add( integral,v.getGlobalIndex(),
								u.getGlobalIndex());
					}
				}
			}
		});
	}
	default void loopRhsViaCell(BiFunction<CT,ST,Double> integralEvaluation,
	                            MatrixFESpace<CT,FT,ET, ST,?,?,?> space){
		List<List<CT>> smallerList = Lists.partition(space.getCells(),space.getCells().size()/PerformanceArguments.getInstance().threadNumber+1);
		Stream<List<CT>> stream = smallerList.stream();
		if(PerformanceArguments.getInstance().parallelizeThreads)
			stream = stream.parallel();
		stream.forEach(smallList->
		{
			for (CT K : smallList)
			{
				for (ST v : space.getShapeFunctionsWithSupportOnCell(K))
				{
					if(space.getFixedNodeIndices().contains(v.getGlobalIndex()))
						continue;
					double integral = integralEvaluation.apply(K,v);
					if(integral != 0)
						space.getRhs().add( integral,v.getGlobalIndex());
				}
			}
		});
	}
	default void loopMatrixViaFace(Function3<FT,ST,ST,Double> integralEvaluation,
	                               MatrixFESpace<CT,FT,ET,ST,?,?,?> space){
		List<List<FT>> smallerList = Lists.partition(space.getFaces(),space.getFaces().size()/PerformanceArguments.getInstance().threadNumber+1);
		Stream<List<FT>> stream = smallerList.stream();
		if(PerformanceArguments.getInstance().parallelizeThreads)
			stream = stream.parallel();
		stream.forEach(smallList->
		{
			for (FT F : smallList)
			{
				System.out.println((int)(50.0*increaseFaceCounter()/space.getFaces().size())+"%");
				for (ST v : space.getShapeFunctionsWithSupportOnFace(F))
				{
					if(space.getFixedNodeIndices().contains(v.getGlobalIndex()))
						continue;
					for (ST u : space.getShapeFunctionsWithSupportOnFace(F))
					{
						if(space.getFixedNodeIndices().contains(u.getGlobalIndex()))
							continue;
						double integral = integralEvaluation.apply(F,u,v);
						if(integral != 0)
							space.getSystemMatrix().add( integral,v.getGlobalIndex(),
								u.getGlobalIndex());
					}
				}
			}
		});
	}
	default void loopRhsViaFace(BiFunction<FT,ST,Double> integralEvaluation,
	                            MatrixFESpace<CT,FT,ET,ST,?,?,?> space){
		List<List<FT>> smallerList = Lists.partition(space.getFaces(),space.getFaces().size()/PerformanceArguments.getInstance().threadNumber+1);
		Stream<List<FT>> stream = smallerList.stream();
		if(PerformanceArguments.getInstance().parallelizeThreads)
			stream = stream.parallel();
		stream.forEach(smallList->
		{
			for (FT F : smallList)
			{
				System.out.println((int)(50.0*increaseFaceCounter()/space.getFaces().size())+"%");
				for (ST v : space.getShapeFunctionsWithSupportOnFace(F))
				{
					if(space.getFixedNodeIndices().contains(v.getGlobalIndex()))
						continue;
					double integral = integralEvaluation.apply(F,v);
					if(integral != 0)
						space.getRhs().add( integral,v.getGlobalIndex());
				}
			}
		});
	}
}
