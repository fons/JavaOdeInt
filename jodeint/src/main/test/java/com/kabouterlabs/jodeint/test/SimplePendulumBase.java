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

import org.bridj.Pointer;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

import static org.bridj.Pointer.pointerToDoubles;

/**
 * Created by fons on 8/30/16.
 */
public class SimplePendulumBase {
    private double t0 = 0.0;
    private double tf = 100.0;
    private double dt = 0.1;
    private int neq = 2;
    private String fn = "lsoda-simple-pendulum.csv";

    private Pointer<Double> qq;
    private double[] q = new double[2];

    static private String filePath(String fn){
        String p = "../data/simple-pendulum/";
        Path path = Paths.get(p);
        if (Files.notExists(path)) {
            return "/.";
        }
        return p + File.separator + fn;

    }

    SimplePendulumBase() {
        q[0] = Math.PI * 999 / 1000.0;
        q[1] = 0.;
        qq = pointerToDoubles(q);
    }
    double t0() {
        return t0;
    }
    double tf() {return tf;}
    double dt() {return dt;}
    int neq() {return neq;}
    Pointer<Double> qq() { return qq;}

    String baseLineResult() {
        return filePath(fn);
    }

}


