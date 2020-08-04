import basic.*;
import com.google.common.base.Stopwatch;
import com.google.common.primitives.Ints;
import linalg.*;
import tensorproduct.*;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
public class Laplace
{
        public static void main(String[] args)
        {


                long startTime = System.nanoTime();

                System.out.println("output start");
                CoordinateVector start = CoordinateVector.fromValues(0,0);
                CoordinateVector end = CoordinateVector.fromValues(1,1);
                int polynomialDegree = 1;
                TPFESpace grid = new TPFESpace(start,end,
                        Ints.asList(20,20),polynomialDegree);
                TPCellIntegral<TPShapeFunction> gg = new TPCellIntegral<>(ScalarFunction.constantFunction(0.01),
                        TPCellIntegral.GRAD_GRAD,
                        false);
                TPFaceIntegral<TPShapeFunction> jj = new TPFaceIntegral<>(ScalarFunction.constantFunction(1000.0),
                        TPFaceIntegral.VALUE_JUMP_VALUE_JUMP, false);
                ArrayList<CellIntegral<TPCell<TPShapeFunction>,TPFace<TPShapeFunction>,TPShapeFunction>> cellIntegrals =
                        new ArrayList<>();
                cellIntegrals.add(gg);
                ArrayList<FaceIntegral<TPCell<TPShapeFunction>,TPFace<TPShapeFunction>,TPShapeFunction>> faceIntegrals = new ArrayList<>();
                faceIntegrals.add(jj);
                TPRightHandSideIntegral<TPShapeFunction> rightHandSideIntegral =
                        new TPRightHandSideIntegral<>(ScalarFunction.constantFunction(1),TPRightHandSideIntegral.VALUE,
                                true);
                ArrayList<RightHandSideIntegral<TPCell<TPShapeFunction>,TPFace<TPShapeFunction>,TPShapeFunction>> rightHandSideIntegrals = new ArrayList<>();
                rightHandSideIntegrals.add(rightHandSideIntegral);
                ArrayList<BoundaryRightHandSideIntegral<TPCell<TPShapeFunction>,TPFace<TPShapeFunction>,TPShapeFunction>> boundaryFaceIntegrals = new ArrayList<>();
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
                IterativeSolver<SparseMatrix> i = new IterativeSolver<>();
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
                Map<CoordinateVector, Double> vals = solut.valuesInPoints(grid.generatePlotPoints(50));
                new PlotFrame(List.of(vals),start,end);
        }
}
