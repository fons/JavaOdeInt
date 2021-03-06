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

import com.kabouterlabs.jodeint.codepack.CodepackLibrary;
import org.bridj.Pointer;
import org.junit.Test;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

import static org.bridj.Pointer.*;
import static org.junit.Assert.assertEquals;


/**
 * Created by fons on 8/5/16.
 */
public class LsodaFullPendulumTest extends SimplePendulumBase
{
    @Test
    public void LsodaFullPendulum() {


        double t  = t0();



        int    lrw = 22 + neq() * Math.max(16, neq()+9);
        int    liw = neq() + 20;
        double rtol = 0.0;
        double atol = 1.0e-12;
        int size = (int) ((tf() - t0()) / dt()) + 2;

        int[]    iwork = new int[liw];
        double[] rwork = new double[lrw];
        double[] q = new double[2];
        q[0] = Math.PI * 999 / 1000.0;
        q[1] = 0.;

        //double[] p = new double[size * (neq + 1)];
        //int stack_size = size * (neq + 1);

        CodepackLibrary.dlsoda_f_callback f = new CodepackLibrary.dlsoda_f_callback() {

            @Override
            public void apply(Pointer<Integer> neq, Pointer<Double> t_, Pointer<Double> q, Pointer<Double> qdot) {
                double alpha = 1;
                qdot.set(0, q.getDoubleAtIndex(1));
                qdot.set(1, -alpha * Math.sin(q.getDoubleAtIndex(0)));
            }
        };
        Pointer<Integer> itolp   = pointerToInt((int) CodepackLibrary.codepack_itol_e.ALL_SCALAR.value());
        Pointer<Integer> itaskp  = pointerToInt((int) CodepackLibrary.codepack_itask_e.NORMAL.value());
        Pointer<Integer> istatep = pointerToInt((int) CodepackLibrary.codepack_istate_in_e.FIRST_CALL.value());
        Pointer<Integer> ioptp   = pointerToInt((int)CodepackLibrary.codepack_iopt_e.NO_OPTIONAL_INPUTS.value());
        Pointer<Integer> jtp     = pointerToInt((int)CodepackLibrary.codepack_jac_type_e.INTERNAL.value());
        Pointer<Double>  atolp   = Pointer.pointerToDouble(atol);
        Pointer<Double>  rtolp   = Pointer.pointerToDouble(rtol);
        Pointer<Integer> neqp    = pointerToInt(neq());
        Pointer<Integer> lrwp    = pointerToInt(lrw);
        Pointer<Integer> liwp    = pointerToInt(liw);
        Pointer<Double> tp       = pointerToDouble(t);

        Pointer<Double> rworkp   = pointerToDoubles(rwork);
        Pointer<Integer> iworkp  = pointerToInts(iwork);
        Pointer<Double> qq       = pointerToDoubles(q);
        Pointer<CodepackLibrary.dlsoda_f_callback> f_func = org.bridj.Pointer.getPointer(f);
        BufferedWriter out;
        String name = new Object(){}.getClass().getEnclosingMethod().getName();

        try
        {
            FileWriter fstream = new FileWriter( PrintStack.path(name + ".csv"), true); //true tells to append data.
            out = new BufferedWriter(fstream);
            while ( tp.get() < tf()) {
                Pointer<Double> tnextp = pointerToDouble(tp.get()+dt());
                CodepackLibrary.dlsoda(f_func,neqp,qq,tp,tnextp, itolp, rtolp, atolp, itaskp, istatep, ioptp, rworkp, lrwp, iworkp, liwp, null, jtp);
                assertEquals("exit code error : " + name, !(istatep.get() < 0),   true);
                out.write(tp.get().toString() + "," + qq.get(0).toString() + "," + qq.get(1) + "\n");
            }
            out.close();
        }
        catch (IOException e)
        {
            System.err.println("Error: " + e.getMessage());
        }
    }

}
