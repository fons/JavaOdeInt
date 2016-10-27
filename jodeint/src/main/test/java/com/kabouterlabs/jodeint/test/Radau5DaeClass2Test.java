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

import com.kabouterlabs.jodeint.cradau5.Cradau5Library;
import org.bridj.Pointer;
import org.junit.Test;

import static org.junit.Assert.assertEquals;

/**
 * Created by fons on 10/24/16.
 */
public class Radau5DaeClass2Test {
    @Test
    public void Radau5DaeClass2() {
        double from = 0.0;
        double to = 10.0;
        double by = 0.1;
        Pointer<Integer> r = Pointer.pointerToInt(2);

        Cradau5Library.radau5_fcn_callback fcn_ = new Cradau5Library.radau5_fcn_callback() {
            @Override
            public void apply(Pointer<Integer> r, Pointer<Double> t, Pointer<Double> y, Pointer<Double> f,
                              Pointer<Double> rpar, Pointer<Integer> ipar) {

                f.set(0, y.get(1));
                f.set(1, y.get(0) - java.lang.Math.cos(t.get()));
            }
        };
        Pointer<Cradau5Library.radau5_fcn_callback> fcn = Pointer.getPointer(fcn_);
        Pointer<Double> t0 = Pointer.pointerToDouble(from);
        Pointer<Double> tend = Pointer.pointerToDouble(from);
        Pointer<Double> y0 = Pointer.allocateDoubles(2);
        y0.set(0, java.lang.Math.cos(0));
        y0.set(1, -java.lang.Math.sin(0));
        Pointer<Double> h = Pointer.pointerToDouble(0.00);
        Pointer<Double> atol = Pointer.pointerToDouble(0.000000001);
        Pointer<Double> rtol = Pointer.pointerToDouble(0.0000000001);
        Pointer<Integer> itol = Pointer.pointerToInt((int) Cradau5Library.radau5_itol_e.ALL_SCALAR.value);

        Pointer<Integer> ijac = Pointer.pointerToInt((int) Cradau5Library.radau5_jacobian_e.INTERNAL.value);

        Cradau5Library.radau5_mas_callback m_ = new Cradau5Library.radau5_mas_callback() {
            @Override

            public void apply(Pointer<Integer> m, Pointer<Double> mas, Pointer<Integer> ldmas, Pointer<Double> rpar, Pointer<Integer> ipar) {

                for (int i = 0; i < m.get()*ldmas.get(); i++) {
                    mas.set(i,0.0);
                }
                mas.set(0, 1.0);
                mas.set(3, 0.0);
            }
        };
        Pointer<Cradau5Library.radau5_mas_callback> mas = Pointer.getPointer(m_);
        Pointer<Integer> imas = Pointer.pointerToInt((int) Cradau5Library.radau5_mass_matrix_e.MASS_USER_PROVIDED.value);
        Pointer<Integer> iout = Pointer.pointerToInt((int) Cradau5Library.radau5_iout_e.NEVER_CALLED.value);
        int ljac = r.get();
        int lmas = r.get();
        int le   = r.get();
        int lw   = r.get()*(ljac + lmas + 3 * le + 12) + 20;
        Pointer<Integer> lwork = Pointer.pointerToInt(lw);
        Pointer<Integer> liwork = Pointer.pointerToInt(r.get()*3+20);
        Pointer<Double> work = Pointer.allocateDoubles(lwork.get());
        Pointer<Integer> iwork = Pointer.allocateInts(liwork.get());
        iwork.set(4,1);//DAE
        iwork.set(5,1);
        iwork.set(6,0);
        Pointer<Integer> idid = Pointer.pointerToInt(0);
        double avg0 = 0.0;
        double avg1 = 0.0;
        int count = 0;
        while (tend.get() < to) {
            count++;
            t0.set(tend.get());
            tend.set(t0.get() + by);
            Cradau5Library.radau5(r, fcn, t0, y0, tend, h, rtol, atol, itol, null, ijac, r, r, mas, imas, r, r, null, iout, work,
                    lwork, iwork, liwork, null, null, idid);
            Double err0 = java.lang.Math.abs(y0.get(0) - java.lang.Math.cos(tend.get()));
            Double err1 = java.lang.Math.abs(y0.get(1) + java.lang.Math.sin(tend.get()));
            avg0 += err0;
            avg1 += err1;
            assertEquals("not a good exit code",(long) idid.get(), (long)1);
            //System.err.print(idid.get().toString() + "   ");
            //System.err.println("result : " + t0.get().toString()+ "," + tend.get().toString() + "," + err0.toString() + "," + y0.get(0).toString()
            //        + "," + err1.toString() + "," + y0.get(1).toString());
        }
        avg0 = avg0 / count;
        avg1 = avg1 / count;
        assertEquals("0 values out of range /s/b cos", 0, avg0, 0.000001);
        assertEquals("1 values out of range /s/b sin", 0, avg1, 0.000001);
    }
}

