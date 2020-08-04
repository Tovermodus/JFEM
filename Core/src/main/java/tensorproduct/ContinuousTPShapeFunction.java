package tensorproduct;

import basic.LagrangeNodeFunctional;
import basic.NodeFunctional;
import basic.ScalarShapeFunction;
import linalg.CoordinateComparator;
import linalg.CoordinateVector;
import linalg.DenseVector;
import linalg.Vector;

import java.util.*;

public class ContinuousTPShapeFunction extends ScalarShapeFunction<TPCell<ContinuousTPShapeFunction>, TPFace<ContinuousTPShapeFunction>,ContinuousTPShapeFunction> implements Comparable<ContinuousTPShapeFunction> {
    
    private Map<TPCell<ContinuousTPShapeFunction>, List<LagrangeBasisFunction1D>> cells;
    private Set<TPFace<ContinuousTPShapeFunction>> faces;
    private final LagrangeNodeFunctional nodeFunctional;
    private final int polynomialDegree;
    private int localIndex;
    public ContinuousTPShapeFunction(TPCell<ContinuousTPShapeFunction> supportCell, int localIndex, int polynomialDegree)
    {
        cells = new TreeMap<>();
        faces = new TreeSet<>();
        this.polynomialDegree = polynomialDegree;
        this.localIndex = localIndex;
        List<LagrangeBasisFunction1D> supportCellFunctions = generateBasisFunctionOnCell(supportCell,
                 localIndex);
        CoordinateVector functionalPoint =
                CoordinateVector.fromValues(supportCellFunctions.stream().mapToDouble(LagrangeBasisFunction1D::getDegreeOfFreedom).toArray());
        nodeFunctional = new LagrangeNodeFunctional(functionalPoint);
        cells.put(supportCell, supportCellFunctions);
        supportCell.addShapeFunction(this);
        checkIfPointOnFace(functionalPoint,supportCell);
//        for(TPFace<ContinuousTPShapeFunction> face: supportCell.faces)
//        {
//            faces.add(face);
//            if (face.isOnFace(functionalPoint))
//            {
//                for (TPCell<ContinuousTPShapeFunction> cellOfFace : face.getCells())
//                {
//                    cells.put(cellOfFace, generateBasisFunctionOnCell(cellOfFace,  functionalPoint));
//                    cellOfFace.addShapeFunction(this);
//                }
//                face.addShapeFunction(this);
//            }
//        }
    }
    private void checkIfPointOnFace(CoordinateVector functionalPoint, TPCell<ContinuousTPShapeFunction> cell)
    {
    
        for(TPFace<ContinuousTPShapeFunction> face: cell.faces)
        {
            if(faces.add(face))
            {
                if (face.isOnFace(functionalPoint))
                {
                    for (TPCell<ContinuousTPShapeFunction> cellOfFace : face.getCells())
                    {
                        cells.put(cellOfFace, generateBasisFunctionOnCell(cellOfFace, functionalPoint));
                        cellOfFace.addShapeFunction(this);
                        checkIfPointOnFace(functionalPoint,cellOfFace);
                    }
                    face.addShapeFunction(this);
                }
            }
        }
    }
    private List<LagrangeBasisFunction1D> generateBasisFunctionOnCell(TPCell<ContinuousTPShapeFunction> cell,
                                                                       int localIndex)
    {
        int[] decomposedLocalIndex = decomposeIndex(cell.getDimension(), polynomialDegree, localIndex);
        List<LagrangeBasisFunction1D> function1Ds = new ArrayList<>();
        for (int i = 0; i < decomposedLocalIndex.length; i++)
        {
            function1Ds.add(new LagrangeBasisFunction1D(polynomialDegree, decomposedLocalIndex[i],
                    cell.cell1Ds.get(i)));
        }
        return function1Ds;
    }
    private List<LagrangeBasisFunction1D> generateBasisFunctionOnCell(TPCell<ContinuousTPShapeFunction> cell,
                                             CoordinateVector functionalPoint)
    {
        List<LagrangeBasisFunction1D> function1Ds = new ArrayList<>();
        for (int i = 0; i < functionalPoint.getLength(); i++)
        {
            function1Ds.add(new LagrangeBasisFunction1D(polynomialDegree, functionalPoint.at(i),
                    cell.cell1Ds.get(i)));
        }
        return function1Ds;
    }
    private static int[] decomposeIndex(int dimension, int polynomialDegree, int localIndex)
    {
        int[] ret = new int[dimension];
        for (int i = 0; i < dimension; i++)
        {
            ret[i] = localIndex % (polynomialDegree+1);
            localIndex = localIndex/(polynomialDegree+1);
        }
        return ret;
    }
    @Override
    public Set<TPCell<ContinuousTPShapeFunction>> getCells() {
        return cells.keySet();
    }

    @Override
    public Set<TPFace<ContinuousTPShapeFunction>> getFaces() {
        return faces;
    }

    @Override
    public NodeFunctional getNodeFunctional() {
        return nodeFunctional;
    }

    @Override
    public void addFace(TPFace<ContinuousTPShapeFunction> face) {
    
        if(faces.add(face))
            face.addShapeFunction(this);
    }

    @Override
    public void addCell(TPCell<ContinuousTPShapeFunction> cell) {
        if(!cells.containsKey(cell))
        {
            cells.put(cell, generateBasisFunctionOnCell(cell, nodeFunctional.getPoint()));
            cell.addShapeFunction(this);
        }
    }


    @Override
    public boolean hasFastEvaluation() {
        return true;
    }
    
    @Override
    public Double value(CoordinateVector pos)
    {
        for(TPCell<ContinuousTPShapeFunction>  c:cells.keySet())
        {
            if(c.isInCell(pos))
                return valueInCell(pos,c);
        }
        return 0.;
    }
    
    @Override
    public Vector gradient(CoordinateVector pos)
    {
        for(TPCell<ContinuousTPShapeFunction>  c:cells.keySet())
        {
            if(c.isInCell(pos))
                return gradientInCell(pos,c);
        }
        return new CoordinateVector(pos.getLength());
    }
    @Override
    public double fastValue(CoordinateVector pos) {
        for(TPCell<ContinuousTPShapeFunction>  c:cells.keySet())
        {
            if(c.isInCell(pos))
                return fastValueInCell(pos,c);
        }
        return 0.;
    }

    @Override
    public double[] fastGradient(CoordinateVector pos) {
        for(TPCell<ContinuousTPShapeFunction>  c:cells.keySet())
        {
            if(c.isInCell(pos))
                return fastGradientInCell(pos,c);
        }
        return new double[pos.getLength()];
    }

    @Override
    public Double valueInCell(CoordinateVector pos, TPCell<ContinuousTPShapeFunction> cell) {
        return fastValueInCell(pos,cell);
    }
    
    @Override
    public Vector gradientInCell(CoordinateVector pos, TPCell<ContinuousTPShapeFunction> cell) {
        return CoordinateVector.fromValues(fastGradientInCell(pos,cell));
    }
    
    @Override
    public double fastValueInCell(CoordinateVector pos, TPCell<ContinuousTPShapeFunction> cell)
    {
        double ret = 1;
        if(cell == null)
            return ret;
        List<? extends Function1D> function1Ds;
        if(cells.containsKey(cell))
        {
            function1Ds = cells.get(cell);
            for (int i = 0; i < pos.getLength(); i++)
            {
                ret *= function1Ds.get(i).value(pos.at(i));
            }
            return ret;
        }
        return 0.;
    }
    
    @Override
    public double[] fastGradientInCell(CoordinateVector pos, TPCell<ContinuousTPShapeFunction> cell)
    {
        double[] ret = new double[pos.getLength()];
        if(cell == null)
            return ret;
        List<? extends Function1D> function1Ds;
        if(cells.containsKey(cell))
        {
            function1Ds = cells.get(cell);
            for (int i = 0; i < pos.getLength(); i++)
            {
                double component = 1;
                for (int j = 0; j < pos.getLength(); j++)
                {
                    if (i == j)
                        component *= function1Ds.get(j).derivative(pos.at(j));
                    else
                        component *= function1Ds.get(j).value(pos.at(j));
                }
                ret[i] = component;
            }
        }
        return ret;
    }
    
    @Override
    public Map<Integer, Double> prolongate(Set<ContinuousTPShapeFunction> refinedFunctions) {
        throw new UnsupportedOperationException();
    }
    
    @Override
    public int compareTo(ContinuousTPShapeFunction o) {
        return CoordinateComparator.comp(nodeFunctional.getPoint(), o.nodeFunctional.getPoint());
    }

    @Override
    public String toString() {
        return "Cell: ".concat(", Node point: ").concat(nodeFunctional.getPoint().toString()).concat(", global Index: ").concat(getGlobalIndex()+"");
    }
    
    @Override
    public boolean equals(Object obj)
    {
        if(obj instanceof ContinuousTPShapeFunction)
            return CoordinateComparator.comp(nodeFunctional.getPoint(),
                    ((ContinuousTPShapeFunction) obj).nodeFunctional.getPoint()) == 0;
        return false;
    }
    
}
