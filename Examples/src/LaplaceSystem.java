import basic.*;
import com.google.common.base.Stopwatch;
import com.google.common.primitives.Ints;
import linalg.CoordinateVector;
import linalg.IterativeSolver;
import linalg.Vector;
import tensorproduct.*;
import tensorproduct.geometry.TPCell;
import tensorproduct.geometry.TPFace;

import java.util.ArrayList;

public class LaplaceSystem
{
        public static void main(String[] args)
        {
        
                PerformanceArguments.PerformanceArgumentBuilder builder =
                        new PerformanceArguments.PerformanceArgumentBuilder();
                builder.build();
                long startTime = System.nanoTime();

                System.out.println("output start");
                CoordinateVector start = CoordinateVector.fromValues(-1,-1);
                CoordinateVector end = CoordinateVector.fromValues(1,1);
                int polynomialDegree = 2;
                TPFESpace grid = new TPFESpace(start,end,
                        Ints.asList(10,5));
                TPCellIntegral<TPShapeFunction> gg = new TPCellIntegral<>(ScalarFunction.constantFunction(1),
                        TPCellIntegral.GRAD_GRAD);
                double penalty = 200000;
                TPFaceIntegral<TPShapeFunction> jj = new TPFaceIntegralViaReferenceFace<>(penalty,
                        TPFaceIntegral.VALUE_JUMP_VALUE_JUMP);
                ArrayList<CellIntegral<TPCell,TPShapeFunction>> cellIntegrals =
                        new ArrayList<>();
                cellIntegrals.add(gg);
                ArrayList<FaceIntegral<TPFace,TPShapeFunction>> faceIntegrals = new ArrayList<>();
                faceIntegrals.add(jj);
                TPRightHandSideIntegral<TPShapeFunction> rightHandSideIntegral =
                        new TPRightHandSideIntegral<>(ScalarFunction.constantFunction(0),TPRightHandSideIntegral.VALUE);
                ArrayList<RightHandSideIntegral<TPCell,TPShapeFunction>> rightHandSideIntegrals = new ArrayList<>();
                rightHandSideIntegrals.add(rightHandSideIntegral);
                TPBoundaryFaceIntegral<TPShapeFunction> bound = new TPBoundaryFaceIntegral<>(new ScalarFunction()
                {
                        @Override
                        public int getDomainDimension()
                        {
                                return 2;
                        }
        
                        @Override
                        public Double value(CoordinateVector pos)
                        {
                                if (Math.abs(pos.x()) == 1||Math.abs(pos.y()) == 1)
                                        return +penalty*2*(1+pos.y())/((3+pos.x())*(3+pos.x())+(1+pos.y())*(1+pos.y()));
                                return (double) 0;
                        }
                },TPBoundaryFaceIntegral.VALUE);
                
                ArrayList<BoundaryRightHandSideIntegral<TPFace,TPShapeFunction>> boundaryFaceIntegrals = new ArrayList<>();
                boundaryFaceIntegrals.add(bound);
                grid.assembleCells();
                grid.assembleFunctions(polynomialDegree);
                grid.initializeSystemMatrix();
                grid.initializeRhs();
                System.out.println("Cell Integrals");
                grid.evaluateCellIntegrals(cellIntegrals,rightHandSideIntegrals);
                System.out.println("Face Integrals");
                grid.evaluateFaceIntegrals(faceIntegrals,boundaryFaceIntegrals);
                System.out.println(((1.0*System.nanoTime() - startTime)/1e9));
                System.out.println("solve system: "+grid.getSystemMatrix().getRows()+"×"+grid.getSystemMatrix().getCols());
                //grid.A.makeParallelReady(12);
                if(grid.getRhs().getLength() < 50)
                {
                        System.out.println(grid.getSystemMatrix());
                        System.out.println(grid.getRhs());
                }
                IterativeSolver i = new IterativeSolver();
                System.out.println("start stopwatch");
                Stopwatch s = Stopwatch.createStarted();
                Vector solution1 = i.solveCG(grid.getSystemMatrix(),grid.getRhs(),1e-3);
                System.out.println(s.elapsed());
                //Vector solution = ((DenseMatrix)grid.getSystemMatrix()).solve(grid.getRhs());
                System.out.println("solved");
                System.out.println(((1.0*System.nanoTime() - startTime)/1e9));
                //grid.A.print_formatted();
                //grid.rhs.print_formatted();
                ScalarFESpaceFunction<TPShapeFunction> solut =
                        new ScalarFESpaceFunction<>(
                                grid.getShapeFunctions(), solution1);
                ScalarFunction referenceSolution = new ScalarFunction()
                {
                        @Override
                        public int getDomainDimension()
                        {
                                return 2;
                        }
        
                        @Override
                        public Double value(CoordinateVector pos)
                        {
                                return 2*(1+pos.y())/((3+pos.x())*(3+pos.x())+(1+pos.y())*(1+pos.y()));
                         }
                };
                PlotWindow p = new PlotWindow();
                p.addPlot(new ScalarPlot2D(solut, grid.generatePlotPoints(20),20));
                p.addPlot(new ScalarPlot2D(solut, grid.generatePlotPoints(50),50));
                //Map<CoordinateVector, Double> vals = solut.valuesInPoints(grid.generatePlotPoints(50));
                //Map<CoordinateVector, Double> refvals = referenceSolution.valuesInPoints(grid.generatePlotPoints(50));
                //new PlotFrame(List.of(vals, refvals),start,end);
        }
}
