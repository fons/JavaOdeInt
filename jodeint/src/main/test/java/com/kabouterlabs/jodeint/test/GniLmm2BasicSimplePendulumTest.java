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

import com.kabouterlabs.jodeint.cgnicodes.CgnicodesLibrary;
import org.bridj.Pointer;
import org.junit.Test;

import java.io.IOException;
import java.io.PrintWriter;

import static org.bridj.Pointer.pointerToDoubles;
import static org.junit.Assert.assertEquals;

/**
 * Created by fons on 8/29/16.
 */
public class GniLmm2BasicSimplePendulumTest extends SimplePendulumBase {
    @Test
    public void GniLmm2BasicSimplePendulum() {

        double t0 = 0.0;
        double tf = 100.0;
        double dt = 0.1;
        int    neq = 1;

        double[] q = new double[1];
        q[0] = Math.PI * 999 / 1000.0;
        double[] qdot = new double[1];
        qdot[0] = 0.;


        CgnicodesLibrary.gnicodes_ode_func f = new CgnicodesLibrary.gnicodes_ode_func(){

            @Override
            public void apply(Pointer<Integer> neq, Pointer<Double> x, Pointer<Double> q, Pointer<Double> f, Pointer<Double> rpar, Pointer<Integer> ipar) {
                double alpha = 1;
                f.set(0, -alpha * Math.sin(q.getDoubleAtIndex(0)));

            }
        };

        Pointer<Double> stack = PrintStack.create(t0,tf,dt,neq * 2);
        Pointer<Double> qq    = pointerToDoubles(q);
        Pointer<Double> qqdot = pointerToDoubles(qdot);
        Pointer<CgnicodesLibrary.gnicodes_ode_func> f_func = org.bridj.Pointer.getPointer(f);


        CgnicodesLibrary.gni_lmm2_basic(stack, qqdot, qq, f_func, neq, t0, tf, dt,CgnicodesLibrary.gni_lmm2_method_e.LMM2_METH_803);
        String name = new Object(){}.getClass().getEnclosingMethod().getName();
        PrintStack.print(stack,t0,tf,dt, 2 * neq,name + ".csv");
        assertEquals("results don't match within expected range : " + name, CsvResultCompare.dataEqual(0.01, stack,t0,tf,dt,2*neq, baseLineResult()),true);
    }

}
