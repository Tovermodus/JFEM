import basic.*;
import com.google.common.base.Stopwatch;
import com.google.common.primitives.Ints;
import linalg.*;
import tensorproduct.*;

import javax.swing.*;
import java.awt.*;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.WindowEvent;
import java.awt.event.WindowListener;
import java.util.ArrayList;
import java.util.Map;
import java.util.OptionalDouble;

public class Laplace
{
        public static void main(String[] args)
        {


                long startTime = System.nanoTime();

                System.out.println("output start");
                CoordinateVector start = CoordinateVector.fromValues(0,0);
                CoordinateVector end = CoordinateVector.fromValues(1,1);
                int polynomialDegree = 2;
                TPFESpace grid = new TPFESpace(start,end,
                        Ints.asList(40,40),polynomialDegree);
                TPCellIntegral gg = new TPCellIntegral(ScalarFunction.constantFunction(1.),TPCellIntegral.GRAD_GRAD,
                        false);
                TPFaceIntegral jj = new TPFaceIntegral(ScalarFunction.constantFunction(1000.0),
                        TPFaceIntegral.VALUE_JUMP_VALUE_JUMP, false);
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
                grid.assembleFunctions(polynomialDegree);
                grid.initializeSystemMatrix();
                grid.initializeRhs();
                grid.evaluateCellIntegrals(cellIntegrals,rightHandSideIntegrals);
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
                        new ScalarFESpaceFunction<TPShapeFunction>(
                                grid.getShapeFunctions(), solution1);
                Map<CoordinateVector, Double> vals = solut.valuesInPoints(grid.generatePlotPoints(50));
                Frame j = new Frame("plot");
        
        
                OptionalDouble maxxo = vals.keySet().stream().mapToDouble(c->c.at(0)).max();
                OptionalDouble minxo = vals.keySet().stream().mapToDouble(c->c.at(0)).min();
                OptionalDouble maxyo = vals.keySet().stream().mapToDouble(c->c.at(1)).max();
                OptionalDouble minyo = vals.keySet().stream().mapToDouble(c->c.at(1)).min();
                double maxx = maxxo.isPresent()? maxxo.getAsDouble(): 0;
                double maxy = maxyo.isPresent()? maxyo.getAsDouble(): 0;
                double minx = minxo.isPresent()? minxo.getAsDouble(): 0;
                double miny = minyo.isPresent()? minyo.getAsDouble(): 0;
                j.setBounds(200,200,
                        (int)(500*(maxx-minx)),(int)(500*(maxy-miny)));
                j.setVisible(true);
                OptionalDouble maxopt = vals.values().stream().mapToDouble(Double::doubleValue).max();
                double max = 0;
                if(maxopt.isPresent())
                        max = maxopt.getAsDouble();
                OptionalDouble minopt = vals.values().stream().mapToDouble(Double::doubleValue).min();
                double min = 0;
                if(minopt.isPresent())
                        min = minopt.getAsDouble();
                j.setTitle("max: "+max+" min: "+min);
                int n = (int)(Math.sqrt(vals.size()));
                //solut.va(100,"/home/tovermodus/plot0.dat");
                System.out.println(((1.0*System.nanoTime() - startTime)/1e9));
                double finalMin = min;
                double finalMax = max;
                j.addMouseListener(new MouseListener()
                                   {
                                           @Override
                                           public void mouseClicked(MouseEvent e)
                                           {
        
                                                   Graphics g = j.getGraphics();
                                                   vals.keySet().stream().forEach(coordinateVector -> {
                                                                           // .getHeight()-100)/n+1);
                                                           double v = (vals.get(coordinateVector) - finalMin)/(finalMax - finalMin);
                                                           g.setColor(new Color((int)(255*v),0,255-(int)(255*v)));
                                                           Vector relc = coordinateVector.sub(start);
                                                           int posx = 50+(int)((j.getWidth()-100)*relc.at(0)/(end.at(0) - start.at(0)));
                                                           int posy = 50+(int)((j.getHeight()-100)*relc.at(1)/(end.at(1) - start.at(1)));
                                                           g.fillRect(posx,posy,(j.getWidth()-100)/n+2,(j.getHeight()-100)/n+2);
                                                   });
                                                   g.setColor(Color.BLACK);
                                                   g.drawString("max: "+ finalMax +" min: "+ finalMin,40,40);
                                                   j.setVisible(true);
                                           }
        
                                           @Override
                                           public void mousePressed(MouseEvent e)
                                           {
                
                                           }
        
                                           @Override
                                           public void mouseReleased(MouseEvent e)
                                           {
                
                                           }
        
                                           @Override
                                           public void mouseEntered(MouseEvent e)
                                           {
                
                                           }
        
                                           @Override
                                           public void mouseExited(MouseEvent e)
                                           {
                                                   System.exit(0);
                                           }
                                   }

                );
        }
}
