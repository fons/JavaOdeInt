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
        String p = "./target/generated-test-data/";
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
