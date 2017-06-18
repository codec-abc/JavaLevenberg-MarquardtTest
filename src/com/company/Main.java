package com.company;

import org.ejml.data.DMatrixRMaj;
import static org.ejml.dense.row.CommonOps_DDRM.*;

public class Main {

    public static void main(String[] args)
    {
        Main main = new Main();
        main.Run();
	// write your code here
    }

    public void Run()
    {
        System.out.println("Hello world");
        myFunc func = new myFunc();
        LevenbergMarquardt lm = new LevenbergMarquardt(func);
        int nbOfPairPoints = 30;
        DMatrixRMaj L /* which is x code*/ = new DMatrixRMaj(nbOfPairPoints, 9 );
        DMatrixRMaj initParam = new DMatrixRMaj(9, 1);
        DMatrixRMaj output = new DMatrixRMaj(nbOfPairPoints ,1);
        lm.optimize(initParam, L, output);
    }

    class myFunc implements LevenbergMarquardt.Function
    {

        @Override
        public void compute
        (
            DMatrixRMaj param, /* x in the paper */
            DMatrixRMaj x, /* L which is fixed */
            DMatrixRMaj y
        )
        {
            //throw new UnsupportedOperationException();
            mult(x, param, y);

        }
    }
}
