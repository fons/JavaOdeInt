package com.kabouterlabs.jodeint.test;

import com.kabouterlabs.jodeint.cradau5.Cradau5Library;
import com.kabouterlabs.jodeint.cradau5.Cradau5Library;
import org.bridj.IntValuedEnum;
import org.bridj.Pointer;
import org.junit.Test;

import java.io.IOException;
import java.io.PrintWriter;

import static org.bridj.Pointer.pointerToDoubles;
import static org.junit.Assert.assertEquals;

/**
 * Created by fons on 8/29/16.
 */
public class Radau5BasicSimplePendulumTest extends SimplePendulumBase {
    @Test
    public void Rdau5BasicSimplePendulum() {


        Cradau5Library.radau5_ode_func f = new Cradau5Library.radau5_ode_func(){

            @Override
            public void apply(Pointer<Integer> neq, Pointer<Double> t_, Pointer<Double> q, Pointer<Double> qdot, Pointer<Double> rpar, Pointer<Integer> ipar) {
                double alpha = 1;
                qdot.set(0, q.getDoubleAtIndex(1));
                qdot.set(1, -alpha * Math.sin(q.getDoubleAtIndex(0)));

            }
        };
        String name = new Object(){}.getClass().getEnclosingMethod().getName();
        Pointer<Double> stack = PrintStack.create(t0(), tf(), dt(), neq());
        Pointer<Cradau5Library.radau5_ode_func> f_func = org.bridj.Pointer.getPointer(f);
        IntValuedEnum<Cradau5Library.radau5_idid_e> error =  Cradau5Library.radau5_basic(stack, qq(), f_func, neq(), t0(), tf(), dt());
        assertEquals("error on :" + name, error,   Cradau5Library.radau5_idid_e.SUCCESS);
        PrintStack.print(stack,t0(),tf(),dt(),neq(),name + ".csv");
        assertEquals("results don't match within expected range : " + name, CsvResultCompare.dataEqual(0.00001, stack,t0(),tf(),dt(),neq(), baseLineResult()),true);
    }

}
