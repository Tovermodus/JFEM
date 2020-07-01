import basic.*;
import com.google.common.base.Stopwatch;
import com.google.common.collect.Iterables;
import linalg.DoubleTensor;
import multigrid.GaussSeidelSmoother;
import multigrid.FESpaceHierarchy;
import scalarspacial.Grid;
import tensorproduct.TPCellIntegral;
import tensorproduct.TPFaceIntegral;
import tensorproduct.TPFESpace;
import tensorproduct.TPRightHandSideIntegral;

import java.util.ArrayList;

public class LaplaceMG
{
        public static void main(String[] args)
        {
                long startTime = System.nanoTime();
                Stopwatch stop = Stopwatch.createStarted();
                System.out.println("output start");

                TPFESpace startGrid = new TPFESpace(0.,0.,1.,1.,3,3,4);
                FESpaceHierarchy fESpaceHierarchy = new FESpaceHierarchy(startGrid);
                CellIntegral gg = new TPCellIntegral(ScalarFunction.constantFunction(1.),TPCellIntegral.GRAD_GRAD);
                TPFaceIntegral jj = new TPFaceIntegral(ScalarFunction.constantFunction(100.0),
                        TPFaceIntegral.VALUE_JUMP_VALUE_JUMP);
                ArrayList<CellIntegral> cellIntegrals = new ArrayList<>();
                cellIntegrals.add(gg);
                ArrayList<FaceIntegral> faceIntegrals = new ArrayList<>();
                faceIntegrals.add(jj);
                TPRightHandSideIntegral rightHandSideIntegral =
                        new TPRightHandSideIntegral(ScalarFunction.oneFunction());
                ArrayList<RightHandSideIntegral> rightHandSideIntegrals = new ArrayList<>();
                rightHandSideIntegrals.add(rightHandSideIntegral);
                ArrayList<BoundaryFaceIntegral> boundaryFaceIntegrals = new ArrayList<>();
                fESpaceHierarchy.addGloballyRefinedLevel();
                fESpaceHierarchy.addGloballyRefinedLevel();
                fESpaceHierarchy.evaluateCellIntegrals(cellIntegrals,rightHandSideIntegrals);
                fESpaceHierarchy.evaluateFaceIntegrals(faceIntegrals,boundaryFaceIntegrals);
                System.out.println(stop.elapsed().toMillis());
                System.out.println("jkjh");
                FESpace g = Iterables.getLast(fESpaceHierarchy.getfESpaces());
                /*
                GaussSeidelSmoother gaussSeidelSmoother = new GaussSeidelSmoother();
                DoubleTensor sol = new DoubleTensor(g.rhs.size());
                for(int i = 0; i < 100; i++)
                {
                        sol = gaussSeidelSmoother.smooth(g,sol,g.rhs);
                        sol.print_formatted();
                }
                FESpaceFunction s = new FESpaceFunction(g.shapeFunctions,sol);
                s.plot(100,"/home/tovermodus/plot0.dat");

                 */


//                for(Face face:gridHierarchy.grids.get(lev).faces)
//                        ((TPFace)face).print();
                DoubleTensor sol = fESpaceHierarchy.multiGridSolve(1e-6,
                100, GaussSeidelSmoother.class,new String[0]);
                //gridHierarchy.grids.get(lev).A.print_formatted(120,105,128,128);
                FESpaceFunction s = new FESpaceFunction(Iterables.getLast(fESpaceHierarchy.getfESpaces()).getShapeFunctions(),sol);
                s.plot(100,"/home/tovermodus/plot0.dat");


        }
}
