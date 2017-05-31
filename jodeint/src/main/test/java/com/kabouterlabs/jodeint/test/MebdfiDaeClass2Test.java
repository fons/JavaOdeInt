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

import com.kabouterlabs.jodeint.cmebdfi.CmebdfiLibrary;
import org.bridj.Pointer;
import org.junit.Test;

import static org.junit.Assert.assertEquals;

/**
 * Created by fons on 10/25/16.
 */
public class MebdfiDaeClass2Test {
    @Test
    public void mebdfiDaeClass2() {
        double from = 0.0;
        double to = 0.4;
        double by = 0.1;
        Pointer<Integer> n = Pointer.pointerToInt(2);

        CmebdfiLibrary.mebdfi_resid_callback fcn_ = new CmebdfiLibrary.mebdfi_resid_callback() {
            @Override
            public void apply(Pointer<Integer> n, Pointer<Double> t, Pointer<Double> y, Pointer<Double> delta, Pointer<Double> yprime,
                              Pointer<Integer> ipar, Pointer<Double> rpar, Pointer<Integer> ier) {
                delta.set(0, yprime.get(0) - y.get(1));
                delta.set(1, y.get(0) - java.lang.Math.cos(t.get()));
            }
        };

        Pointer<CmebdfiLibrary.mebdfi_resid_callback> resid = Pointer.getPointer(fcn_);

        CmebdfiLibrary.mebdfi_pderv_callback p_ = new CmebdfiLibrary.mebdfi_pderv_callback() {
            @Override
            public void apply(Pointer<Double> doublePtr1, Pointer<Double> doublePtr2, Pointer<Double> doublePtr3, Pointer<Integer> intPtr1,
                              Pointer<Double> doublePtr4, Pointer<Integer> intPtr2, Pointer<Double> doublePtr5, Pointer<Integer> intPtr3,
                              Pointer<Double> doublePtr6, Pointer<Integer> intPtr4) {

            }
        };

        Pointer<CmebdfiLibrary.mebdfi_pderv_callback> pderv = Pointer.getPointer(p_);
        Pointer<Double> t0 = Pointer.pointerToDouble(from);
        Pointer<Double> tout = Pointer.pointerToDouble(from);
        Pointer<Double> tend = Pointer.pointerToDouble(to);

        Pointer<Double> y0 = Pointer.allocateDoubles(2);
        y0.set(0, java.lang.Math.cos(0));
        y0.set(1, -java.lang.Math.sin(0));

        Pointer<Double> yprime = Pointer.allocateDoubles(2);
        yprime.set(0, -java.lang.Math.sin(0));
        yprime.set(1, -java.lang.Math.cos(0));

        Pointer<Integer> mf = Pointer.pointerToInt((int)CmebdfiLibrary.mebdfi_method_e.BDF_INTERNAL_FULL_JAC.value);
        Pointer<Integer> idid = Pointer.pointerToInt((int)CmebdfiLibrary.mebdfi_idid_in_e.FIRST_CALL.value);
        Pointer<Integer> lout  = Pointer.pointerToInt(6);
        Pointer<Integer> maxder   = Pointer.pointerToInt(7);

        Pointer<Double> h0 = Pointer.pointerToDouble(0.00001);
        Pointer<Double> atol = Pointer.pointerToDouble(0.00001);
        Pointer<Double> rtol = Pointer.pointerToDouble(0.00001);
        Pointer<Integer> itol = Pointer.pointerToInt((int) CmebdfiLibrary.mebdfi_itol_e.ALL_SCALAR.value);


        Pointer<Integer> ipar = Pointer.allocateInt();
        Pointer<Double>  rpar = Pointer.allocateDouble();

        Pointer<Integer> mbnd = Pointer.allocateInts(4);
        mbnd.set(0, n.get());
        mbnd.set(1, n.get());
        mbnd.set(2, n.get());
        mbnd.set(3, n.get());

        Pointer<Integer> lwork = Pointer.pointerToInt(4 + 32 * n.get() + 2* n.get() * mbnd.get(3));
        Pointer<Integer> liwork = Pointer.pointerToInt(n.get() + 14);
        Pointer<Double>  work = Pointer.allocateDoubles(lwork.get());
        Pointer<Integer> iwork = Pointer.allocateInts(liwork.get());

        iwork.set(0,1);//DAE
        iwork.set(1,1);
        iwork.set(2,0);
        iwork.set(13,5000);

        Pointer<Integer> ier = Pointer.pointerToInt(0);
        double avg0 = 0.0;
        double avg1 = 0.0;
        int count = 0;

        int index = 0;
        int max   = java.lang.Math.abs((int)((to - from)/by));
        while (index < max) {
            //t0.set(from + index * by);
            t0.set(tout.get());
            tout.set(t0.get()+by);
            index++;
            count++;

            ///System.err.println(t0.get().toString() + "," + tout.get().toString());
            CmebdfiLibrary.mebdfi(n, t0, h0, y0, yprime, tout, tend, mf, idid, lout, lwork, work, liwork, iwork,mbnd, maxder, itol, rtol,atol, rpar,
                    ipar, pderv, resid, ier);

            if (idid.get() == 1) {
                idid.set(0);
                CmebdfiLibrary.mebdfi(n, t0, h0, y0, yprime, tout, tend, mf, idid, lout, lwork, work, liwork, iwork,mbnd, maxder, itol, rtol,atol, rpar,
                        ipar, pderv, resid, ier);
            }
            if (idid.get() == -6) {
                System.err.println("===> " + iwork.get(13).toString());
            }
            //System.err.println("===> " + iwork.get(13).toString());
            //System.err.println("===> " + work.get(1).toString());
            assertEquals("not a good exit code",(long) 0, (long) idid.get());
            Double err0 = java.lang.Math.abs(y0.get(0) - java.lang.Math.cos(tout.get()));
            Double err1 = java.lang.Math.abs(y0.get(1) + java.lang.Math.sin(tout.get()));
            avg0 += err0;
            avg1 += err1;
            //System.err.print(idid.get().toString() + "   ");
            //System.err.println("result : " + t0.get().toString()+ "," + tout.get().toString() + "," + err0.toString() + "," + y0.get(0).toString()
            //        + "," + err1.toString() + "," + y0.get(1).toString());
            idid.set(1);
            //work.set(1,0.0001);
        }
        avg0 = avg0 / count;
        avg1 = avg1 / count;
        assertEquals("0 values out of range /s/b cos", 0, avg0, 0.000001);
        assertEquals("1 values out of range /s/b sin", 0, avg1, 0.0001);
    }
}
