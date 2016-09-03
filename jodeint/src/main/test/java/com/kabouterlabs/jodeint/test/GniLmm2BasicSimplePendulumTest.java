package com.kabouterlabs.jodeint.test;

import com.kabouterlabs.jodeint.cgnicodes.CgnicodesLibrary;
import org.bridj.Pointer;
import org.junit.Test;

import java.io.IOException;
import java.io.PrintWriter;

import static org.bridj.Pointer.pointerToDoubles;
import static org.junit.Assert.assertEquals;

/**
 * Created by fons on 8/29/16.
 */
public class GniLmm2BasicSimplePendulumTest extends SimplePendulumBase {
    @Test
    public void GniLmm2BasicSimplePendulum() {

        double t0 = 0.0;
        double tf = 100.0;
        double dt = 0.1;
        int    neq = 1;

        double[] q = new double[1];
        q[0] = Math.PI * 999 / 1000.0;
        double[] qdot = new double[1];
        qdot[0] = 0.;


        CgnicodesLibrary.gnicodes_ode_func f = new CgnicodesLibrary.gnicodes_ode_func(){

            @Override
            public void apply(Pointer<Integer> neq, Pointer<Double> x, Pointer<Double> q, Pointer<Double> f, Pointer<Double> rpar, Pointer<Integer> ipar) {
                double alpha = 1;
                f.set(0, -alpha * Math.sin(q.getDoubleAtIndex(0)));

            }
        };

        Pointer<Double> stack = PrintStack.create(t0,tf,dt,neq * 2);
        Pointer<Double> qq    = pointerToDoubles(q);
        Pointer<Double> qqdot = pointerToDoubles(qdot);
        Pointer<CgnicodesLibrary.gnicodes_ode_func> f_func = org.bridj.Pointer.getPointer(f);


        CgnicodesLibrary.gni_lmm2_basic(stack, qqdot, qq, f_func, neq, t0, tf, dt,CgnicodesLibrary.gni_lmm2_method_e.LMM2_METH_803);
        String name = new Object(){}.getClass().getEnclosingMethod().getName();
        PrintStack.print(stack,t0,tf,dt, 2 * neq,name + ".csv");
        assertEquals("results don't match within expected range : " + name, CsvResultCompare.dataEqual(0.00001, stack,t0,tf,dt,2*neq, baseLineResult()),true);
    }

}
