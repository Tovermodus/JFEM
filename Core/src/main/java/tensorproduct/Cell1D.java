package tensorproduct;

import org.jetbrains.annotations.NotNull;

import java.util.ArrayList;

public class Cell1D implements Comparable<Cell1D>
{
        private final double start;
        private final double end;
        public ArrayList<LagrangeBasisFunction1D> shapefunctions;
        public Cell1D(double start, double end)
        {
                this.start = start;
                this.end = end;
        }
        public Cell1D(Cell1D cell)
        {
                this(cell.start, cell.end);
        }
        
        public double[][] distributeQuadrature(QuadratureRule1D quadratureRule)
        {
                double[][] pointsWeights = new double[2][quadratureRule.length()];
                for(int i = 0; i < quadratureRule.length(); i++)
                {
                        pointsWeights[0][i] = positionOnGrid(quadratureRule.getReferencePoints().get(i));
                        pointsWeights[1][i] = quadratureRule.getReferenceWeights().get(i)*length();
                }
                return pointsWeights;
        }
        public double getStart()
        {
                return start;
        }
        
        public double getEnd()
        {
                return end;
        }
        
        public double length()
        {
                return (end - start);
        }

        public double center()
        {
                return (end + start) / 2;
        }

        public boolean isInCell(double pos)                //in [0,1]
        {
                return (pos >= start && pos <= end);
        }

        public double positionOnReferenceCell(double pos)
        {
                return (pos - start) / length();
        }

        public double positionOnGrid(double pospp)
        {
                return pospp * length() + start;
        }

        public double jacobiDeterminant(double pos)
        {
                return 1.0 / length();
        }

        public void distributeFunctions(int polynomialDegree)
        {
                for(int  i = 0; i < polynomialDegree + 1; i++)
                {
                        shapefunctions.add(new LagrangeBasisFunction1D(polynomialDegree, i, this));
                }
        }
        
        @Override
        public String toString()
        {
                return "Cell: start" + start + ", end: "+ end;
        }
        
        public void print()
        {
                System.out.println(this);
        }
        
        @Override
        public int compareTo(@NotNull Cell1D o)
        {
                if(getStart() < o.getStart())
                        return -1;
                if(getStart() > o.getStart())
                        return 1;
                return 0;
        }
}
