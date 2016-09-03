package com.kabouterlabs.jodeint.test;



import org.bridj.Pointer;
import org.bridj.PointerIO;
import org.junit.Test;
import com.kabouterlabs.jodeint.cquadpack.CquadpackLibrary;
import org.bridj.IntValuedEnum;

import static org.junit.Assert.assertEquals;

/**
 * Created by fons on 8/7/16.
 */
public class QagseBasicTest {
    @Test
    public void QagseBasic1() {
        CquadpackLibrary.cquadpack_ode_func f = new CquadpackLibrary.cquadpack_ode_func() {
            @Override
            public double apply(Pointer<Double> x) {
                return 4.0 * x.get() * x.get() - 5.0 * x.get();
            }
        };
        Pointer<CquadpackLibrary.cquadpack_ode_func> f_ptr = Pointer.getPointer(f);

        double start  = 0.0;
        double end    = 10.0;
        double expected = 4000/3.0 - 250.0;
        Pointer<Double> result_ptr = Pointer.allocate(PointerIO.getDoubleInstance());

        IntValuedEnum<CquadpackLibrary.quadpack_errno_e> errno = CquadpackLibrary.qagse_basic(f_ptr, start, end, result_ptr);
        assertEquals("error on the odeint ", errno,  CquadpackLibrary.quadpack_errno_e.SUCCESS);
        assertEquals("result of integration",expected, result_ptr.get().doubleValue(), 0.0000001);
        //System.err.println(" result : " + result_ptr.get().toString() + "  expected " + expected + " errno :" + errno.toString());
    }

    @Test
    public void QagseBasic2() {
        CquadpackLibrary.cquadpack_ode_func f = new CquadpackLibrary.cquadpack_ode_func() {
            @Override
            public double apply(Pointer<Double> x) {
                return Math.log(x.get()) /Math.sqrt(x.get());
            }
        };
        Pointer<CquadpackLibrary.cquadpack_ode_func> f_ptr = Pointer.getPointer(f);

        double start  = 0.0;
        double end    = 1.0;
        double expected = -4.0;
        Pointer<Double> result_ptr = Pointer.allocate(PointerIO.getDoubleInstance());

        IntValuedEnum<CquadpackLibrary.quadpack_errno_e> errno = CquadpackLibrary.qagse_basic(f_ptr, start, end, result_ptr);
        assertEquals("error on the odeint ", errno,  CquadpackLibrary.quadpack_errno_e.SUCCESS);
        assertEquals("result of integration", expected, result_ptr.get().doubleValue(), 0.0000001);

        //System.err.println(" result : " + result_ptr.get().toString() + "  expected " + expected + " errno :" + errno.toString());
    }

}
