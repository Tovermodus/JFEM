package basic;

import io.vavr.Tuple2;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.function.BiConsumer;
import java.util.function.BinaryOperator;
import java.util.function.Function;
import java.util.function.Supplier;
import java.util.stream.Collector;

public class MapCollector<S, T>
	implements Collector<Tuple2<S, T>, Map<S, T>, Map<S, T>>
{
	@Override
	public Supplier<Map<S, T>> supplier()
	{
		return HashMap::new;
	}
	
	@Override
	public BiConsumer<Map<S, T>, Tuple2<S, T>> accumulator()
	{
		return (m, t) -> m.put(t._1, t._2);
	}
	
	@Override
	public BinaryOperator<Map<S, T>> combiner()
	{
		return (stMap, stMap2) ->
		{
			stMap.putAll(stMap2);
			return stMap;
		};
	}
	
	@Override
	public Function<Map<S, T>, Map<S, T>> finisher()
	{
		return stMap -> stMap;
	}
	
	@Override
	public Set<Characteristics> characteristics()
	{
		return new HashSet<>();
	}
}
