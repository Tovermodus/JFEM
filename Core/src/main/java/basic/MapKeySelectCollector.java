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

public class MapKeySelectCollector<Kprev, Kaft, T>
	implements Collector<Tuple2<Kprev, T>, Map<Kaft, T>, Map<Kaft, T>>
{
	
	final java.util.function.Function<Kprev, Kaft> keySelector;
	
	public MapKeySelectCollector(final java.util.function.Function<Kprev, Kaft> keySelector)
	{
		this.keySelector = keySelector;
	}
	
	@Override
	public Supplier<Map<Kaft, T>> supplier()
	{
		return HashMap::new;
	}
	
	@Override
	public BiConsumer<Map<Kaft, T>, Tuple2<Kprev, T>> accumulator()
	{
		return (m, t) -> m.put(keySelector.apply(t._1), t._2);
	}
	
	@Override
	public BinaryOperator<Map<Kaft, T>> combiner()
	{
		return (stMap, stMap2) ->
		{
			stMap.putAll(stMap2);
			return stMap;
		};
	}
	
	@Override
	public Function<Map<Kaft, T>, Map<Kaft, T>> finisher()
	{
		return stMap -> stMap;
	}
	
	@Override
	public Set<Characteristics> characteristics()
	{
		return new HashSet<>();
	}
}
