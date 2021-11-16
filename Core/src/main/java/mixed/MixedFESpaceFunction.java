package mixed;

import basic.Cell;
import basic.Face;
import basic.ScalarFunction;
import basic.VectorFunction;
import linalg.CoordinateDenseMatrix;
import linalg.CoordinateMatrix;
import linalg.CoordinateVector;
import linalg.Vector;

import java.util.HashMap;
import java.util.Map;
import java.util.Optional;
import java.util.concurrent.ConcurrentHashMap;

public class MixedFESpaceFunction<MF extends MixedShapeFunction<CT, FT, ?, ?>, CT extends Cell<CT, FT>,
	FT extends Face<CT, FT>>
	extends MixedFunctionOnCells<CT, FT>
{
	final double diam;
	private final HashMap<MF, Double> coefficients;
	Map<MF, Double> pressureCoefficients;
	Map<MF, Double> velocityCoefficients;
	
	public MixedFESpaceFunction(final MF[] functions, final double[] coefficients)
	{
		super();
		diam = functions[0]
			.getCells()
			.stream()
			.findAny()
			.orElseThrow(() -> new IllegalArgumentException("Function Has No Cells"))
			.diam();
		assert (functions.length == coefficients.length);
		this.coefficients = new HashMap<>();
		for (int i = 0; i < functions.length; i++)
		{
			this.coefficients.put(functions[i], coefficients[i]);
		}
		initializeFunctionSets();
	}
	
	public MixedFESpaceFunction(final Map<Integer, MF> functions, final Vector coefficients)
	{
		assert (functions.size() == coefficients.size());
		diam = functions
			.values()
			.stream()
			.findAny()
			.orElseThrow(() -> new IllegalArgumentException("Has no Functions"))
			.getCells()
			.stream()
			.findAny()
			.orElseThrow(() -> new IllegalArgumentException("Function Has No Cells"))
			.diam();
		;
		this.coefficients = new HashMap<>();
		for (final Map.Entry<Integer, MF> function : functions.entrySet())
		{
			this.coefficients.put(function.getValue(), coefficients.at(function.getKey()));
		}
		initializeFunctionSets();
	}
	
	private void initializeFunctionSets()
	{
		pressureCoefficients = new ConcurrentHashMap<>();
		velocityCoefficients = new ConcurrentHashMap<>();
		for (final MF shapeFunction : coefficients.keySet())
		{
			if (shapeFunction.hasPressureFunction())
				pressureCoefficients.put(shapeFunction, coefficients.get(shapeFunction));
			else velocityCoefficients.put(shapeFunction, coefficients.get(shapeFunction));
		}
	}
	
	@Override
	public int getDomainDimension()
	{
		return coefficients.keySet()
		                   .iterator()
		                   .next()
		                   .getDomainDimension();
	}
	
	@Override
	public ScalarFunction getPressureFunction()
	{
		final MixedFESpaceFunction<MF, CT, FT> me = this;
		return new ScalarFunction()
		{
			@Override
			public int getDomainDimension()
			{
				return me.getDomainDimension();
			}
			
			@Override
			public Double value(final CoordinateVector pos)
			{
				return me.value(pos)
				         .getPressure();
			}
			
			@Override
			public CoordinateVector gradient(final CoordinateVector pos)
			{
				return me.gradient(pos)
				         .getPressureGradient();
			}
		};
	}
	
	@Override
	public VectorFunction getVelocityFunction()
	{
		final MixedFESpaceFunction<MF, CT, FT> me = this;
		return new VectorFunction()
		{
			@Override
			public int getRangeDimension()
			{
				return me.value(new CoordinateVector(getDomainDimension()))
				         .getLength();
			}
			
			@Override
			public int getDomainDimension()
			{
				return me.getDomainDimension();
			}
			
			@Override
			public CoordinateVector value(final CoordinateVector pos)
			{
				return me.value(pos)
				         .getVelocity();
			}
			
			@Override
			public CoordinateMatrix gradient(final CoordinateVector pos)
			{
				return me.gradient(pos)
				         .getVelocityGradient();
			}
		};
	}
	
	@Override
	public boolean hasPressureFunction()
	{
		return true;
	}
	
	@Override
	public boolean hasVelocityFunction()
	{
		return true;
	}
	
	@Override
	public MixedValue value(final CoordinateVector pos)
	{
		final PressureValue pf = new PressureValue(pressureCoefficients
			                                           .entrySet()
			                                           .stream()
			                                           .parallel()
			                                           .filter(entry -> entry
				                                           .getKey()
				                                           .hasPressureFunction())
			                                           .mapToDouble(entry -> entry
				                                           .getKey()
				                                           .value(pos)
				                                           .getPressure() * entry.getValue())
			                                           .sum());
		final VelocityValue vf = velocityCoefficients
			.entrySet()
			.stream()
			.parallel()
			.filter(entry -> entry.getKey()
			                      .hasVelocityFunction())
			.map(entry -> (VelocityValue) (entry.getKey()
			                                    .value(pos)).mul(entry.getValue()))
			.reduce(new VelocityValue(getDomainDimension()), (a, b) -> (VelocityValue) (a.add(b)));
		return pf.add(vf);
	}
	
	@Override
	public MixedGradient gradient(final CoordinateVector pos)
	{
		final Optional<CoordinateVector> pressureGradient = pressureCoefficients
			.entrySet()
			.stream()
			.parallel()
			.map(entry -> entry.getKey()
			                   .gradient(pos)
			                   .getPressureGradient()
			                   .mul(entry.getValue()))
			.reduce(CoordinateVector::add);
		PressureGradient pf = new PressureGradient(pos.mul(0));
		if (pressureGradient.isPresent()) pf = new PressureGradient(pressureGradient.get());
		final CoordinateDenseMatrix velocityGradient = velocityCoefficients
			.entrySet()
			.stream()
			.parallel()
			.map(entry -> entry.getKey()
			                   .gradient(pos)
			                   .getVelocityGradient()
			                   .mul(entry.getValue()))
			.reduce(new CoordinateDenseMatrix(getDomainDimension() + 1, getDomainDimension() + 1),
			        CoordinateDenseMatrix::add);
		final VelocityGradient vf = new VelocityGradient(velocityGradient);
		return pf.add(vf);
	}
}
