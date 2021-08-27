package basic;

import linalg.CoordinateDenseMatrix;
import linalg.CoordinateMatrix;
import linalg.CoordinateVector;
import linalg.DenseMatrix;

import java.util.List;

public class FEBaseTransformation<F extends Function<valueT, gradientT, hessianT>, NT extends NodeFunctional<valueT, gradientT, hessianT>, valueT, gradientT, hessianT>
{
	List<F> originalBasis;
	List<NT> nodeFunctionals;
	public DenseMatrix transformationMatrix;
	
	public FEBaseTransformation(final List<F> functions, final List<NT> nodeFunctionals)
	{
		if (nodeFunctionals.size() != functions.size()) throw new IllegalArgumentException(
			"need same number of nodefunctionals as basisfunctions" + functions.size() + " " + nodeFunctionals
				.size());
		originalBasis = functions;
		this.nodeFunctionals = nodeFunctionals;
		final DenseMatrix nodeValuesOfFunctions = new DenseMatrix(nodeFunctionals.size(), functions.size());
		for (int i = 0; i < nodeFunctionals.size(); i++)
			for (int j = 0; j < functions.size(); j++)
				nodeValuesOfFunctions.set(nodeFunctionals.get(i).evaluate(functions.get(j)), j, i);
		transformationMatrix = nodeValuesOfFunctions.inverse();
	}
	
	public double scalarBasisFunctionValue(final NT nodeFunctional, final CoordinateVector pos)
	{
		return scalarBasisFunctionValue(nodeFunctionals.indexOf(nodeFunctional), pos);
	}
	
	public CoordinateVector scalarBasisFunctionGradient(final NT nodeFunctional, final CoordinateVector pos)
	{
		return scalarBasisFunctionGradient(nodeFunctionals.indexOf(nodeFunctional), pos);
	}
	
	public CoordinateVector vectorBasisFunctionValue(final NT nodeFunctional, final CoordinateVector pos)
	{
		return vectorBasisFunctionValue(nodeFunctionals.indexOf(nodeFunctional), pos);
	}
	
	public CoordinateMatrix vectorBasisFunctionGradient(final NT nodeFunctional, final CoordinateVector pos)
	{
		return vectorBasisFunctionGradient(nodeFunctionals.indexOf(nodeFunctional), pos);
	}
	
	public double vectorBasisFunctionDivergence(final NT nodeFunctional, final CoordinateVector pos)
	{
		return vectorBasisFunctionDivergence(nodeFunctionals.indexOf(nodeFunctional), pos);
	}
	
	public double scalarBasisFunctionValue(final int index, final CoordinateVector pos)
	{
		double ret = 0;
		for (int i = 0; i < originalBasis.size(); i++)
		{
			ret += (double) originalBasis.get(i).value(pos) * transformationMatrix.at(index, i);
		}
		return ret;
	}
	
	public CoordinateVector scalarBasisFunctionGradient(final int index, final CoordinateVector pos)
	{
		CoordinateVector ret = new CoordinateVector(pos.getLength());
		for (int i = 0; i < originalBasis.size(); i++)
			ret = ret.add(((CoordinateVector) originalBasis.get(i).gradient(pos)).mul(
				transformationMatrix.at(index, i)));
		return ret;
	}
	
	public CoordinateVector vectorBasisFunctionValue(final int index, final CoordinateVector pos)
	{
		CoordinateVector ret = new CoordinateVector(pos.getLength());
		for (int i = 0; i < originalBasis.size(); i++)
		{
			//System.out.println(ret);
			ret = ret.add(((CoordinateVector) originalBasis.get(i).value(pos)).mul(
				transformationMatrix.at(index, i)));
		}
		return ret;
	}
	
	public CoordinateMatrix vectorBasisFunctionGradient(final int index, final CoordinateVector pos)
	{
		CoordinateMatrix ret = new CoordinateDenseMatrix(pos.getLength(), pos.getLength());
		for (int i = 0; i < originalBasis.size(); i++)
		{
			//System.out.println(originalBasis.get(i).gradient(pos)+" orig");
			ret = ret.add(((CoordinateMatrix) originalBasis.get(i).gradient(pos)).mul(
				transformationMatrix.at(index, i)));
		}
		//System.out.println(ret);
		return ret;
	}
	
	public Double vectorBasisFunctionDivergence(final int index, final CoordinateVector pos)
	{
		double ret = 0;
		for (int i = 0; i < originalBasis.size(); i++)
		{
			//System.out.println(((VectorFunction)originalBasis.get(i)).divergence(pos)+" orig");
			ret += ((VectorFunction) originalBasis.get(i)).divergence(pos) * transformationMatrix.at(index,
			                                                                                         i);
		}
		//System.out.println(ret);
		return ret;
	}
}
