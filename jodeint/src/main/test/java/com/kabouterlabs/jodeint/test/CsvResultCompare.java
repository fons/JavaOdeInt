
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

/**
 * Created by fons on 8/31/16.
 */

import au.com.bytecode.opencsv.CSVReader;
import org.bridj.Pointer;


import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;

import java.util.List;

public class CsvResultCompare
{
    static private double MULT = Math.log(0.50);
    static private double ZEROISH = Math.log10(Double.MIN_VALUE);

    static private List readCsv(String fn)
    {
        CSVReader reader;
        try {
            reader = new CSVReader(new FileReader(fn));
        } catch (FileNotFoundException e) {
            e.printStackTrace();
            return null;
        }
        try {
            return reader.readAll();
        } catch (IOException e) {
            e.printStackTrace();
        }

        return null;
    }

    static private double[] csvToList(List<String[]> csv, int size) {
        double[] cstack = new double[size];
        int index = 0;

        for (String[] line : csv) {
            for (String item : line) {
                if (index < size && item.length() > 0) {
                    cstack[index] = Double.parseDouble(item);
                    index++;
                }
            }
        }

        return cstack;
    }

    static double diff(double l, double r)
    {
        if ( l == r) return 0;
        //double n = Math.log(Math.abs(l - r));
        //double d = Math.log(Math.abs(l + r));
        //double diff = Math.exp(n - d + MULT);
        double diff = Math.abs(l-r);
        return diff;
    }
    static public Boolean dataEqual(double tolerance, Pointer<Double> stackptr, double t0, double tf, double dt, int neq, String fn)
    {

        List<String[]> csv = readCsv(fn);
        int stack_size = (int) ((neq + 1) * ((tf - t0) / dt)) + 2;
        double[] cstack = csvToList(csv, stack_size);
        double[] stack = (double[]) stackptr.getArray();
        double lavg = 0;
        double avg  = 0;
        int    ak = 0;
        for (int i = 0; i < stack_size; i++) {

            double ard = diff(cstack[i], stack[i]);
            //System.err.println(((Integer)i).toString() + "," + ((Double)cstack[i]).toString() + "," + ((Double)stack[i]).toString() +"," + ((Double)ard).toString());
            double lard = 0;//ZEROISH;
            if ( ard != 0.0) {
                lard = Math.log10(ard);
                lavg = ((ak * lavg) + lard)/ (ak + 1.0);
                ak++;
            }
            avg = ((i * avg) + ard)/ (i + 1.0);

        }
        return (lavg < tolerance && Math.log10(avg) < tolerance);
    }
}
