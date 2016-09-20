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


import com.kabouterlabs.jodeint.crkf45.Crkf45Library;
import org.bridj.IntValuedEnum;
import org.bridj.Pointer;
import org.junit.Test;

import static org.junit.Assert.assertEquals;

/**
 * Created by fons on 8/29/16.
 */
public class Rkf45BasicSimplePendulumTest extends SimplePendulumBase {
    @Test
    public void Rkf45BasicSimplePendulum() {

        Crkf45Library.rkf45_ode_func f = new Crkf45Library.rkf45_ode_func(){
            @Override
            public void apply(Pointer<Double> t_, Pointer<Double> q, Pointer<Double> qdot) {
                double alpha = 1;
                qdot.set(0, q.getDoubleAtIndex(1));
                qdot.set(1, -alpha * Math.sin(q.getDoubleAtIndex(0)));

            }
        };
        String name = new Object(){}.getClass().getEnclosingMethod().getName();
        Pointer<Double> stack = PrintStack.create(t0(),tf(),dt(),neq());
        Pointer<Crkf45Library.rkf45_ode_func> f_func = org.bridj.Pointer.getPointer(f);

        IntValuedEnum<Crkf45Library.rkf45_retval_e> error =  Crkf45Library.rkf45_basic(stack, qq(), f_func, neq(), t0(), tf(), dt());
        assertEquals("error on " + name , error,   Crkf45Library.rkf45_retval_e.SUCCESS);
        PrintStack.print(stack,t0(),tf(),dt(),neq(),name + ".csv");
        assertEquals("results don't match within expected range : " + name, CsvResultCompare.dataEqual(0.00001, stack,t0(),tf(),dt(),neq(), baseLineResult()),true);
    }

}
