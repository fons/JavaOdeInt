package com.kabouterlabs.jodeint.test;

import com.kabouterlabs.jodeint.cdop853.Cdop853Library;
import org.bridj.IntValuedEnum;
import org.bridj.Pointer;
import org.junit.Test;


import static org.junit.Assert.assertEquals;

/**
 * Created by fons on 8/26/16.
 */
public class Dop853BasicSimplePendulumTest extends SimplePendulumBase {
    @Test
    public void Dop853BasicSimplePendulum() {


        Cdop853Library.dop853_ode_func f = new Cdop853Library.dop853_ode_func(){

            @Override
            public void apply(Pointer<Integer> neq, Pointer<Double> t_, Pointer<Double> q, Pointer<Double> qdot, Pointer<Double> rpar, Pointer<Integer> ipar) {
                double alpha = 1;
                qdot.set(0, q.getDoubleAtIndex(1));
                qdot.set(1, -alpha * Math.sin(q.getDoubleAtIndex(0)));

            }
        };
        String name = new Object(){}.getClass().getEnclosingMethod().getName();
        Pointer<Double> stack = PrintStack.create(t0(),tf(),dt(),neq());
        Pointer<Cdop853Library.dop853_ode_func> f_func = org.bridj.Pointer.getPointer(f);
        IntValuedEnum<Cdop853Library.dop853_idid_e> error =  Cdop853Library.dop853_basic(stack, qq(), f_func, neq(), t0(), tf(), dt());
        assertEquals("error on : " + name, error,   Cdop853Library.dop853_idid_e.SUCCESS);

        PrintStack.print(stack,t0(),tf(),dt(),neq(),name + ".csv");
        assertEquals("results don't match within expected range : " + name, CsvResultCompare.dataEqual(0.00001, stack,t0(),tf(),dt(),neq(), baseLineResult()),true);
    }

}
