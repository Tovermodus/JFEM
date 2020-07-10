import basic.*;
import com.google.common.primitives.Ints;
import linalg.CoordinateVector;
import linalg.DenseMatrix;
import linalg.DenseVector;
import tensorproduct.*;

import java.util.ArrayList;

public class Laplace
{
        public static void main(String[] args)
        {


                long startTime = System.nanoTime();

                System.out.println("output start");
                TPFESpace grid = new TPFESpace(CoordinateVector.fromValues(0,0),CoordinateVector.fromValues(1,1),
                        Ints.asList(2,2),1);
                TPCellIntegral gg = new TPCellIntegral(TPCellIntegral.GRAD_GRAD);
                TPFaceIntegral jj = new TPFaceIntegral(ScalarFunction.constantFunction(1000.0),
                        TPFaceIntegral.VALUE_JUMP_VALUE_JUMP, true);
                ArrayList<CellIntegral<TPCell<TPShapeFunction>,TPFace<TPShapeFunction>,TPShapeFunction>> cellIntegrals =
                        new ArrayList<>();
                cellIntegrals.add(gg);
                ArrayList<FaceIntegral<TPCell<TPShapeFunction>,TPFace<TPShapeFunction>,TPShapeFunction>> faceIntegrals = new ArrayList<>();
                faceIntegrals.add(jj);
                TPRightHandSideIntegral rightHandSideIntegral =
                        new TPRightHandSideIntegral(ScalarFunction.constantFunction(1),TPRightHandSideIntegral.VALUE,
                                true);
                ArrayList<RightHandSideIntegral<TPCell<TPShapeFunction>,TPFace<TPShapeFunction>,TPShapeFunction>> rightHandSideIntegrals = new ArrayList<>();
                rightHandSideIntegrals.add(rightHandSideIntegral);
                ArrayList<BoundaryRightHandSideIntegral<TPCell<TPShapeFunction>,TPFace<TPShapeFunction>,TPShapeFunction>> boundaryFaceIntegrals = new ArrayList<>();
                grid.assembleCells();
                grid.assembleFunctions(1);
                grid.initializeSystemMatrix();
                grid.initializeRhs();
                grid.evaluateCellIntegrals(cellIntegrals,rightHandSideIntegrals);
                grid.evaluateFaceIntegrals(faceIntegrals,boundaryFaceIntegrals);
                System.out.println(((1.0*System.nanoTime() - startTime)/1e9));
                System.out.println("solve system: "+grid.getSystemMatrix().getRows()+"×"+grid.getSystemMatrix().getCols());
                //grid.A.makeParallelReady(12);
                System.out.println(grid.getSystemMatrix());
                System.out.println(grid.getRhs());
                DenseVector solution = ((DenseMatrix)grid.getSystemMatrix()).solve(grid.getRhs());
                System.out.println("solved");
                System.out.println(((1.0*System.nanoTime() - startTime)/1e9));
                //grid.A.print_formatted();
                //grid.rhs.print_formatted();
                ScalarFESpaceFunction<TPShapeFunction> solut =
                        new ScalarFESpaceFunction<TPShapeFunction>(
                                grid.getShapeFunctions(), solution);
                solut.valuesInPoints(grid.generatePlotPoints(0.5));
                //solut.va(100,"/home/tovermodus/plot0.dat");
                System.out.println(((1.0*System.nanoTime() - startTime)/1e9));
        }
}
