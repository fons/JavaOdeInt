package com.kabouterlabs.jodeint.test;


import com.kabouterlabs.jodeint.crkf45.Crkf45Library;
import org.bridj.IntValuedEnum;
import org.bridj.Pointer;
import org.junit.Test;

import static org.junit.Assert.assertEquals;

/**
 * Created by fons on 8/29/16.
 */
public class Rkf45BasicSimplePendulumTest extends SimplePendulumBase {
    @Test
    public void Rkf45BasicSimplePendulum() {

        Crkf45Library.rkf45_ode_func f = new Crkf45Library.rkf45_ode_func(){
            @Override
            public void apply(Pointer<Double> t_, Pointer<Double> q, Pointer<Double> qdot) {
                double alpha = 1;
                qdot.set(0, q.getDoubleAtIndex(1));
                qdot.set(1, -alpha * Math.sin(q.getDoubleAtIndex(0)));

            }
        };
        String name = new Object(){}.getClass().getEnclosingMethod().getName();
        Pointer<Double> stack = PrintStack.create(t0(),tf(),dt(),neq());
        Pointer<Crkf45Library.rkf45_ode_func> f_func = org.bridj.Pointer.getPointer(f);

        IntValuedEnum<Crkf45Library.rkf45_retval_e> error =  Crkf45Library.rkf45_basic(stack, qq(), f_func, neq(), t0(), tf(), dt());
        assertEquals("error on " + name , error,   Crkf45Library.rkf45_retval_e.SUCCESS);
        PrintStack.print(stack,t0(),tf(),dt(),neq(),name + ".csv");
        assertEquals("results don't match within expected range : " + name, CsvResultCompare.dataEqual(0.00001, stack,t0(),tf(),dt(),neq(), baseLineResult()),true);
    }

}
