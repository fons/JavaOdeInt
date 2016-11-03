package com.kabouterlabs.jodeint.test;

import com.kabouterlabs.jodeint.cdaspk.CdaspkLibrary;
import org.bridj.Pointer;
import org.junit.Test;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

/**
 * Created by fons on 10/28/16.
 */
public class DaspkDaeClass2Test {
    @Test
    public void daspkDaeClass2() {
        double from = 0.0;
        double to = 10.0;
        double by = 0.1;


        CdaspkLibrary.daspk_res_callback fcn_ = new CdaspkLibrary.daspk_res_callback() {
            @Override
            public void apply(Pointer<Double> t, Pointer<Double> y, Pointer<Double> yprime, Pointer<Double> cj, Pointer<Double> delta,
                              Pointer<Integer> ires, Pointer<Double> rpar, Pointer<Integer> ipar) {
                delta.set(0, yprime.get(0) - y.get(1));
                delta.set(1, y.get(0) - java.lang.Math.cos(t.get()));

            }


        };
        Pointer<CdaspkLibrary.daspk_res_callback> res = Pointer.getPointer(fcn_);
        Pointer<Integer> neq = Pointer.pointerToInt(2);


        Pointer<Double> t0 = Pointer.pointerToDouble(from);
        Pointer<Double> tout = Pointer.pointerToDouble(from);


        Pointer<Double> y0 = Pointer.allocateDoubles(2);
        y0.set(0, java.lang.Math.cos(0));
        y0.set(1, -java.lang.Math.sin(0));

        Pointer<Double> yprime = Pointer.allocateDoubles(2);
        yprime.set(0, -java.lang.Math.sin(0));
        yprime.set(1, -java.lang.Math.cos(0));

        Pointer<Integer> info = Pointer.allocateInts(40 + neq.get());

        //info.set(20,1);
        //info.set(21,1);
        //info.set(22,0);

        Pointer<Integer> idid = Pointer.pointerToInt(1);

        Pointer<Double> atol = Pointer.pointerToDouble(0.00001);
        Pointer<Double> rtol = Pointer.pointerToDouble(0.00000001);
        Pointer<Integer> lrw = Pointer.pointerToInt(50 + 9 * neq.get() + neq.get() * neq.get());
        Pointer<Integer> liw = Pointer.pointerToInt(neq.get() + 40);

        Pointer<Double>  rwork = Pointer.allocateDoubles(lrw.get());
        Pointer<Integer> iwork = Pointer.allocateInts(liw.get());


        double avg0 = 0.0;
        double avg1 = 0.0;
        int count = 0;
        int max_reps = 10;
        int index = 0;
        int max   = java.lang.Math.abs((int)((to - from)/by));
        while (index < max) {
            int reps = 0;
            boolean continue_rep = true;
            index++;
            count++;
            //t0.set(tout.get());
            tout.set(t0.get()+by);
            //info.set(0,0);
            while (continue_rep) {
                CdaspkLibrary.daspk(res, neq, t0, y0, yprime, tout, info, rtol, atol, idid, rwork, lrw, iwork, liw, null, null, null, null);
                if (idid.get() == -1 && reps < max_reps) {
                    reps++;
                    info.set(0, 1);
                    //System.err.println("===> " + rwork.get(2).toString());
                    //System.err.println(t0.get().toString() + "," + tout.get().toString());
                } else {
                    continue_rep = false;
                }
                if (idid.get() > 0) {
                    //System.err.println("===> " + rwork.get(2).toString());
                    continue_rep = false;
                    Double err0 = java.lang.Math.abs(y0.get(0) - java.lang.Math.cos(tout.get()));
                    Double err1 = java.lang.Math.abs(y0.get(1) + java.lang.Math.sin(tout.get()));
                    //System.err.println("result : " + t0.get().toString()+ "," + tout.get().toString() + "," + err0.toString() + "," + y0.get(0).toString()
                    //        + "," + err1.toString() + "," + y0.get(1).toString());
                }
            }

            assertTrue("not a good exit code :" + idid.get().toString()  ,(long) idid.get() > 0);
            Double err0 = java.lang.Math.abs(y0.get(0) - java.lang.Math.cos(tout.get()));
            Double err1 = java.lang.Math.abs(y0.get(1) + java.lang.Math.sin(tout.get()));
            avg0 += err0;
            avg1 += err1;
            //System.err.print(idid.get().toString() + "   ");
            //System.err.println("result : " + t0.get().toString()+ "," + tout.get().toString() + "," + err0.toString() + "," + y0.get(0).toString()
            //        + "," + err1.toString() + "," + y0.get(1).toString());
            idid.set(1);
        }
        avg0 = avg0 / count;
        avg1 = avg1 / count;
        assertEquals("0 values out of range /s/b cos", 0, avg0, 0.000001);
        assertEquals("1 values out of range /s/b sin", 0, avg1, 0.0001);
    }}
