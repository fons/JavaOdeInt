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


import com.kabouterlabs.jodeint.czvode.CzvodeLibrary;
import com.kabouterlabs.jodeint.czvode.zvode_complex_s;
import org.bridj.IntValuedEnum;
import org.bridj.Pointer;
import org.junit.Test;

import java.io.IOException;
import java.io.PrintWriter;

import static org.bridj.Pointer.pointerToDoubles;
import static org.junit.Assert.assertEquals;

/**
 * Created by fons on 8/5/16.
 */
public class ZvodeBasicSimplePendulumTest extends SimplePendulumBase {
    @Test
    public void  ZvodeBasicSimplePendulum()
    {

        zvode_complex_s[] q = new zvode_complex_s[2];
        q[0] = new zvode_complex_s();
        q[1] = new zvode_complex_s();
        q[0].dr(Math.PI * 999 / 1000.0);
        q[0].di(0.0);
        q[1].dr(0.0);
        q[1].di(0.0);


        CzvodeLibrary.czvode_ode_func ff5 = new CzvodeLibrary.czvode_ode_func() {
            @Override
            public void apply(Pointer<Integer> neq, Pointer<Double> t_, Pointer<zvode_complex_s> y, Pointer<zvode_complex_s> ydot,
                              Pointer<zvode_complex_s> rpar, Pointer<Integer> ipar) {

                double alpha = 1;


                zvode_complex_s v0 = new zvode_complex_s();
                v0.dr(y.get(1).dr());
                v0.di(y.get(1).di());

                ydot.set(0, v0);

                zvode_complex_s v = new zvode_complex_s();
                v.dr(-alpha * Math.sin(y.get(0).dr()));
                v.di(0.0);

                ydot.set(1,v);


            }

        };
        String name = new Object(){}.getClass().getEnclosingMethod().getName();
        Pointer<Double> stack = PrintStack.create(t0(), tf(),dt(), 2 * neq());

        Pointer<zvode_complex_s> qq = Pointer.allocateArray(zvode_complex_s.class, 2);
        qq.set(0, q[0]);
        qq.set(1, q[1]);

        Pointer<CzvodeLibrary.czvode_ode_func> f_func = org.bridj.Pointer.getPointer(ff5);
        IntValuedEnum<CzvodeLibrary.czvode_ode_err_e> error = CzvodeLibrary.zvode_basic(stack, qq, f_func, neq(), t0(), tf(), dt(), CzvodeLibrary.czvode_method_e.ADAMS);
        assertEquals("error on : " + name, error,   CzvodeLibrary.czvode_ode_err_e.SUCCESS);
        PrintStack.print(stack, t0(),tf(),dt(), 2 * neq(),  name + ".csv");

    }

}

