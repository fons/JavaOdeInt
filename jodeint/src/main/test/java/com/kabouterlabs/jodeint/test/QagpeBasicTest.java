/*
 * Copyright (c) 2016.
 * JodeInt Developers
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 */

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

        //////IntValuedEnum<CquadpackLibrary.quadpack_errno_e> errno = CquadpackLibrary.qagpe_basic(f_ptr, start,end,npts2, points_ptr, result_ptr);
        ///assertEquals("error on the qagpe_basic ", CquadpackLibrary.quadpack_errno_e.SUCCESS,errno);
        ///assertEquals("result of integration",expected, result_ptr.get().doubleValue(), 0.0000001);
        //System.err.println(" result : " + result_ptr.get().toString() + "  expected " + expected + " errno :" + errno.toString());
    }
}
