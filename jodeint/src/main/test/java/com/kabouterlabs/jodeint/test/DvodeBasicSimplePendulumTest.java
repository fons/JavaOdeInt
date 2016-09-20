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

import com.kabouterlabs.jodeint.cdvode.CdvodeLibrary;


import org.junit.Test;

import static org.junit.Assert.assertEquals;

import org.bridj.*;

import java.io.IOException;
import java.io.PrintWriter;

import static org.bridj.Pointer.*;

/**
 * Created by fons on 7/27/16.
 */
public class DvodeBasicSimplePendulumTest extends SimplePendulumBase {
    @Test
    public void  DvodeBasicSimplePendulum()
    {


        CdvodeLibrary.cdvode_ode_func ff5 = new CdvodeLibrary.cdvode_ode_func() {

            public void apply(Pointer<Integer> neq, Pointer<Double> t_, Pointer<Double> q, Pointer<Double> qdot) {
                double alpha = 1;
                qdot.set(0, q.getDoubleAtIndex(1));
                qdot.set(1, -alpha * Math.sin(q.getDoubleAtIndex(0)));
            }
        };
        Pointer<Double> stack = PrintStack.create(t0(), tf(), dt(), neq());

        Pointer<CdvodeLibrary.cdvode_ode_func> f_func = org.bridj.Pointer.getPointer(ff5);

        String name = new Object(){}.getClass().getEnclosingMethod().getName();
        IntValuedEnum<CdvodeLibrary.cdvode_ode_err_e> error = CdvodeLibrary.dvode_basic(stack, qq(), f_func, neq(), t0(), tf(), dt(),
                CdvodeLibrary.cdvode_method_e.ADAMS);
        assertEquals("error on " + name, error,   CdvodeLibrary.cdvode_ode_err_e.SUCCESS);


        PrintStack.print(stack,t0(),tf(),dt(),neq(),name + ".csv");

        assertEquals("results don't match within expected range : " + name, CsvResultCompare.dataEqual(0.00001, stack,t0(),tf(),dt(),neq(), baseLineResult()),true);
    }

}
