package com.kabouterlabs.jodeint.test;

import com.kabouterlabs.jodeint.codepack.CodepackLibrary;
import org.bridj.IntValuedEnum;
import org.bridj.Pointer;
import org.junit.Test;

import java.io.IOException;
import java.io.PrintWriter;

import static org.bridj.Pointer.pointerToDoubles;
import static org.junit.Assert.assertEquals;

/**
 * Created by fons on 7/29/16.
 */
public class LsodaBasicSimplePendulumTest extends SimplePendulumBase {
    @Test
    public void LsodaSimplePendulum() {

        String name = new Object(){}.getClass().getEnclosingMethod().getName();

        CodepackLibrary.codepack_ode_func f = new CodepackLibrary.codepack_ode_func() {

            @Override
            public void apply(Pointer<Integer> neq, Pointer<Double> t_, Pointer<Double> q, Pointer<Double> qdot) {
                double alpha = 1;
                qdot.set(0, q.getDoubleAtIndex(1));
                qdot.set(1, -alpha * Math.sin(q.getDoubleAtIndex(0)));
            }
        };

        Pointer<Double> stack = PrintStack.create(t0(),tf(),dt(),neq());
        Pointer<CodepackLibrary.codepack_ode_func> f_func = org.bridj.Pointer.getPointer(f);
        IntValuedEnum<CodepackLibrary.codepack_ode_err_e> error =  CodepackLibrary.lsoda_basic(stack, qq(), f_func, neq(), t0(), tf(), dt());

        assertEquals("exit code error : " + name, error,   CodepackLibrary.codepack_ode_err_e.SUCCESS);

        PrintStack.print(stack,t0(),tf(),dt(),neq(),name + ".csv");

        assertEquals("results don't match within expected range : " + name, CsvResultCompare.dataEqual(0.00001, stack,t0(),tf(),dt(),neq(), baseLineResult()),true);

    }

}
