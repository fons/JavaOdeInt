package com.kabouterlabs.jodeint.test;

import com.kabouterlabs.jodeint.cquadpack.CquadpackLibrary;
import org.bridj.IntValuedEnum;
import org.bridj.Pointer;
import org.bridj.PointerIO;
import org.junit.Test;

import static org.junit.Assert.assertEquals;

/**
 * Created by fons on 8/8/16.
 */
public class QawseBasicTest {
    @Test
    public void QawseBasic() {
        CquadpackLibrary.cquadpack_ode_func f = new CquadpackLibrary.cquadpack_ode_func() {
            @Override
            public double apply(Pointer<Double> x) {
                if (x.get() > 0) {
                    double v1 = Math.log(x.get());
                    double v2 = 1.0 + v1*v1;
                    return 1.0/(v2*v2);
                }
              return 0.0;
            }
        };
        Pointer<CquadpackLibrary.cquadpack_ode_func> f_ptr = Pointer.getPointer(f);

        double start    = 0.0;
        double end      = 1.0;
        double alpha    = 0.0;
        double beta     = 0.0;

        double ci = 0.33740392290096813466;

        double si = 0.94608307036718301494;
        double expected = 0.5* ( ci * Math.sin (1.0) + (Math.PI / 2.0- si ) * Math.cos (1.0) - 1.0);
        Pointer<Double> result_ptr = Pointer.allocate(PointerIO.getDoubleInstance());
        //System.err.println(CquadpackLibrary.log_weight_to_string(CquadpackLibrary.quadpack_log_weight_function_e.LOGW_1));
        IntValuedEnum<CquadpackLibrary.quadpack_errno_e> errno = CquadpackLibrary.qawse_basic(f_ptr, start,end,alpha,beta,
                    CquadpackLibrary.quadpack_log_weight_function_e.LOGW_2, result_ptr);
        assertEquals("error on the odeint ", errno,  CquadpackLibrary.quadpack_errno_e.SUCCESS);
        assertEquals("result of integration",expected, result_ptr.get().doubleValue(), 0.0000001);
        //System.err.println(" result : " + result_ptr.get().toString() + "  expected " + expected + " errno :" + errno.toString());

    }
}
