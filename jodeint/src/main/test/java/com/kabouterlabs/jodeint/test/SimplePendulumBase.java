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


