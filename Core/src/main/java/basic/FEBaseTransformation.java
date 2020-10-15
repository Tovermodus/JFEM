package basic;

import linalg.CoordinateMatrix;
import linalg.CoordinateVector;
import linalg.DenseMatrix;

import java.util.List;

public class FEBaseTransformation<F extends Function<valueT,	gradientT,	hessianT>,
	NT extends NodeFunctional<F,
	valueT,
	gradientT,hessianT>, valueT,
	gradientT, hessianT>
{
	List<F> originalBasis;
	List<NT> nodeFunctionals;
	public DenseMatrix transformationMatrix;
	public  FEBaseTransformation(List<F> functions, List<NT> nodeFunctionals)
	{
		if(nodeFunctionals.size() != functions.size())
			throw new IllegalArgumentException("need same number of nodefunctionals as basisfunctions"+ functions.size()+ " "+ nodeFunctionals.size());
		originalBasis = functions;
		this.nodeFunctionals = nodeFunctionals;
		DenseMatrix nodeValuesOfFunctions = new DenseMatrix(nodeFunctionals.size(), functions.size());
		for(int i = 0; i < nodeFunctionals.size(); i++)
			for(int j = 0; j < functions.size(); j++)
				nodeValuesOfFunctions.set(nodeFunctionals.get(i).evaluate(functions.get(j)),i,j);
		transformationMatrix = nodeValuesOfFunctions.inverse();
	}
	public double scalarBasisFunctionValue(NT nodeFunctional, CoordinateVector pos)
	{
		return scalarBasisFunctionValue(nodeFunctionals.indexOf(nodeFunctional), pos);
	}
	public CoordinateVector scalarBasisFunctionGradient(NT nodeFunctional, CoordinateVector pos)
	{
		return scalarBasisFunctionGradient(nodeFunctionals.indexOf(nodeFunctional), pos);
	}
	public CoordinateVector vectorBasisFunctionValue(NT nodeFunctional, CoordinateVector pos)
	{
		return vectorBasisFunctionValue(nodeFunctionals.indexOf(nodeFunctional), pos);
	}
	public CoordinateMatrix vectorBasisFunctionGradient(NT nodeFunctional, CoordinateVector pos)
	{
		return vectorBasisFunctionGradient(nodeFunctionals.indexOf(nodeFunctional), pos);
	}
	public double scalarBasisFunctionValue(int index, CoordinateVector pos)
	{
		double ret = 0;
		for(int i = 0; i < originalBasis.size(); i++)
			ret += (double)originalBasis.get(i).value(pos)*transformationMatrix.at(index,i);
		return ret;
	}
	public CoordinateVector scalarBasisFunctionGradient(int index, CoordinateVector pos)
	{
		CoordinateVector ret = new CoordinateVector(pos.getLength());
		for(int i = 0; i < originalBasis.size(); i++)
			ret =
				ret.add (((CoordinateVector)originalBasis.get(i).gradient(pos)).mul(transformationMatrix.at(index
					,i)));
		return ret;
	}
	public CoordinateVector vectorBasisFunctionValue(int index, CoordinateVector pos)
	{
		CoordinateVector ret = new CoordinateVector(pos.getLength());
		for(int i = 0; i < originalBasis.size(); i++)
			ret =
				ret.add (((CoordinateVector)originalBasis.get(i).value(pos)).mul(transformationMatrix.at(index
					,i)));
		return ret;
	}
	public CoordinateMatrix vectorBasisFunctionGradient(int index, CoordinateVector pos)
	{
		CoordinateMatrix ret = new CoordinateMatrix(pos.getLength(), pos.getLength());
		for(int i = 0; i < originalBasis.size(); i++)
			ret =
				ret.add (((CoordinateMatrix)originalBasis.get(i).gradient(pos)).mul(transformationMatrix.at(index
					,i)));
		return ret;
	}
}
