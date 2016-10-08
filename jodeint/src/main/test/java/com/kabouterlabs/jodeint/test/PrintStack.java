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
import java.io.IOException;
import java.io.PrintWriter;

import static org.junit.Assert.assertEquals;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.nio.file.Path;
/**
 * Created by fons on 8/30/16.
 */
public class PrintStack {


    static  String path(String fn) {
        String p = "../target/test/generated-test-data/";
        Path path = Paths.get(p);
        if (Files.notExists(path)) {
            try {
                Files.createDirectories(path);
            } catch (IOException e) {
                e.printStackTrace();
                return File.separator + fn;
            }
        }
        return p + File.separator + fn;
    }
    static Pointer<Double> create(double t0, double tf, double dt, int neq) {
        int size = (int) ((tf - t0) / dt) + 2;
        double[] p = new double[size * (neq + 1)];
        return Pointer.pointerToDoubles(p);
    }

    static void print(Pointer<Double> stack, double t0, double tf, double dt, int neq, String fn)
    {
        //System.err.println("saved data to " + fn);
        int size = (int) ((tf - t0) / dt) + 2;
        int stack_size = size * (neq + 1);
        try {
            PrintWriter writer1 = new PrintWriter(path(fn), "UTF-8");
            double[] n = (double[]) stack.getArray();
            for (int i = 0; i < stack_size; i++) {
                writer1.print(n[i]);
                writer1.print(",");
                if ((i + 1) % (neq + 1) == 0) {
                    writer1.println("");
                }
            }

            writer1.close();
        } catch (IOException e) {
            e.printStackTrace();
            assertEquals("io error", 0, 1);
        }

    }
}
