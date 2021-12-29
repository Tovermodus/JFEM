package mixed;

import basic.Cell;
import basic.Face;
import basic.ScalarFunctionOnCells;
import basic.VectorFunctionOnCells;
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
	implements MixedFunctionOnCells<CT, FT>
{
	final double diam;
	private final HashMap<MF, Double> coefficients;
	Map<MF, Double> pressureCoefficients;
	Map<MF, Double> velocityCoefficients;
	final int dimension;
	
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
		dimension = this.coefficients.keySet()
		                             .iterator()
		                             .next()
		                             .getDomainDimension();
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
		dimension = this.coefficients.keySet()
		                             .iterator()
		                             .next()
		                             .getDomainDimension();
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
		return dimension;
	}
	
	@Override
	public ScalarFunctionOnCells<CT, FT> getPressureFunction()
	{
		final MixedFESpaceFunction<MF, CT, FT> me = this;
		return new ScalarFunctionOnCells<CT, FT>()
		{
			@Override
			public Double valueInCell(final CoordinateVector pos, final CT cell)
			{
				return me.valueInCell(pos, cell)
				         .getPressure();
			}
			
			@Override
			public CoordinateVector gradientInCell(final CoordinateVector pos, final CT cell)
			{
				return me.gradientInCell(pos, cell)
				         .getPressureGradient();
			}
			
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
	public VectorFunctionOnCells<CT, FT> getVelocityFunction()
	{
		final MixedFESpaceFunction<MF, CT, FT> me = this;
		return new VectorFunctionOnCells<CT, FT>()
		{
			@Override
			public CoordinateVector valueInCell(final CoordinateVector pos, final CT cell)
			{
				return me.valueInCell(pos, cell)
				         .getVelocity();
			}
			
			@Override
			public CoordinateMatrix gradientInCell(final CoordinateVector pos, final CT cell)
			{
				return me.gradientInCell(pos, cell)
				         .getVelocityGradient();
			}
			
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
		final MixedValue pf = new PressureValue(pressureCoefficients
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
			                                        .sum(), dimension);
		final MixedValue vf = velocityCoefficients
			.entrySet()
			.stream()
			.parallel()
			.filter(entry -> entry.getKey()
			                      .hasVelocityFunction())
			.map(entry -> entry.getKey()
			                   .value(pos)
			                   .mul(entry.getValue()))
			.reduce(new VelocityValue(getDomainDimension()), MixedValue::add);
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
	
	@Override
	public MixedValue valueInCell(final CoordinateVector pos, final CT cell)
	{
		return value(pos);
	}
	
	@Override
	public MixedGradient gradientInCell(final CoordinateVector pos, final CT cell)
	{
		return gradient(pos);
	}
}
