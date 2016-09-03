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
public class QagpeBasicTest {
    @Test
    public void QagpeBasic() {
        CquadpackLibrary.cquadpack_ode_func f = new CquadpackLibrary.cquadpack_ode_func() {
            @Override
            public double apply(Pointer<Double> x) {
                    double v1 = x.get();
                    double v2 = v1*v1;
                    return v1*v2* Math.log(Math.abs((v2-2.0)*(v2-1.0)));
            }
        };
        Pointer<CquadpackLibrary.cquadpack_ode_func> f_ptr = Pointer.getPointer(f);

        double start    = 0.0;
        double end      = 3.0;
        int    npts2    = 4;

        double[] points = new double[npts2];
        points[0] = 1.0;
        points[1] = Math.sqrt(2.0);

        Pointer<Double> points_ptr = Pointer.pointerToDoubles(points);

        double expected = 61.0*Math.log(2)+77.0*Math.log(7.0)/4.0 - 27.0;
        Pointer<Double> result_ptr = Pointer.allocate(PointerIO.getDoubleInstance());

        IntValuedEnum<CquadpackLibrary.quadpack_errno_e> errno = CquadpackLibrary.qagpe_basic(f_ptr, start,end,npts2, points_ptr, result_ptr);
        assertEquals("error on the qagpe_basic ", CquadpackLibrary.quadpack_errno_e.SUCCESS,errno);
        assertEquals("result of integration",expected, result_ptr.get().doubleValue(), 0.0000001);
        //System.err.println(" result : " + result_ptr.get().toString() + "  expected " + expected + " errno :" + errno.toString());
    }
}
