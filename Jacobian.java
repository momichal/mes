import static java.lang.StrictMath.sqrt;

public class Jacobian {
    private double[][][][] J;
    private double [][][][] inverseMatrix;
    private double[][] derKsi;
    private double[][] derEta;
    private Element[] elements;
    private Node[] nodes;

    private double det;
    private int dt;

    private double[] dNdx;
    private double[] dNdy;

    private double[][] Hx;
    private double[][] Hy;
    private double[][] H;
    private double[][] Hsum;
    private double[][] Hlocal;
    private double[][] Hglobal;
    private double[][] Hfinal;

    private double[][] N;
    private double[][] Nsum;

    private double[][] Cglobal;

    private double[][] warBrzeg;
    private double[] P;
    private double[] Pglobal;

    private int EL;
    private int ND;

    private double[] T0;

    public Jacobian(Grid grid){
        EL = grid.getEL().length;
        ND = grid.getND().length;

        J = new double[EL][4][2][2];
        inverseMatrix = new double[EL][4][2][2];
        derKsi = new double[4][4];
        derEta = new double[4][4];
        elements = new Element[EL];
        nodes = new Node[ND];

        dNdx = new double[4];
        dNdy = new double[4];

        Hx = new double[4][4];
        Hy = new double[4][4];
        H = new double[4][4];
        Hsum = new double[4][4];
        Hlocal = new double[4][4];
        Hglobal = new double[ND][ND];
        Hfinal = new double[ND][ND];

        N = new double[4][4];
        Nsum = new double[4][4];

        Cglobal = new double[ND][ND];

        warBrzeg = new double[4][4];
        P = new double[4];
        Pglobal = new double[ND];

        T0 = new double[ND];
    }

    public void jacobian(UniversalElement uniEl, Grid grid) {
        derKsi = uniEl.getDerivativeKsi();
        derEta = uniEl.getDerivativeEta();
        elements = grid.getEL();
        nodes = grid.getND();

        double x, y;
        for(int i = 1; i <= elements.length; i++) {
            System.out.println("Jacobian - element " + i);
            int [] elementID = elements[i-1].getID();
            for(int p = 0; p < 4; p++){
                //System.out.println("Punkt " + p);
                for(int j = 0; j < 4; j++){
                    //System.out.println("Wezel " + j);
                    x = nodes[elementID[j]-1].getX();    /* -1 bo odwołuje się np. do nodes[1] a tak na prawdę to nodes[0] */
                    y = nodes[elementID[j]-1].getY();
                    //System.out.println("x: " + x + ", y: " + y + ", pochKsi: " + derKsi[p][j] + ", pochEta: " + derEta[p][j]);
                    J[i-1][p][0][0] += derKsi[p][j] * x;
                    J[i-1][p][0][1] += derKsi[p][j] * y;
                    J[i-1][p][1][0] += derEta[p][j] * x;
                    J[i-1][p][1][1] += derEta[p][j] * y;
                }
                det = J[i-1][p][0][0] * J[i-1][p][1][1] - (J[i-1][p][0][1] * J[i-1][p][1][0]);

                inverseMatrix[i-1][p][0][0] = 1/det * J[i-1][p][1][1];
                inverseMatrix[i-1][p][0][1] = 1/det * (-J[i-1][p][0][1]);
                inverseMatrix[i-1][p][1][0] = 1/det * (-J[i-1][p][1][0]);
                inverseMatrix[i-1][p][1][1] = 1/det * J[i-1][p][0][0];
                System.out.println();
                showJacobian(i, p, det, inverseMatrix);
            }
            System.out.println();
        }
    }

    public void transformationJacobian(Grid grid){
        double alfa = grid.getGlobalData().getAlfa();
        double B = grid.getGlobalData().getB();
        int nB = grid.getGlobalData().get_nB() - 1; //liczba wezlow na szerokosc -1
        double delta = B/nB; //długosc boku

        elements = grid.getEL();
        nodes = grid.getND();

        double wspK = grid.getGlobalData().getK();

        for(int i = 1; i <= elements.length; i++) {
            System.out.println("Transformation Jacobian - element " + i);
            for (int p = 0; p < 4; p++) {
                System.out.println("the integration point: " + (p + 1));
                for (int j = 0; j < 4; j++) {
                    dNdx[j] = inverseMatrix[i - 1][p][0][0] * derKsi[p][j] + inverseMatrix[i - 1][p][0][1] * derEta[p][j];
                    dNdy[j] = inverseMatrix[i - 1][p][1][0] * derKsi[p][j] + inverseMatrix[i - 1][p][1][1] * derEta[p][j];
                }
                System.out.print("dN/dx : ");
                for (int j = 0; j < 4; j++) {
                    System.out.print(dNdx[j] + "; ");
                }
                System.out.println();
                System.out.print("dN/dy : ");
                for (int j = 0; j < 4; j++) {
                    System.out.print(dNdy[j] + "; ");
                }
                System.out.println();

                for (int k = 0; k < 4; k++) {
                    for (int j = 0; j < 4; j++) {
                        Hx[k][j] = (dNdx[j] * dNdx[k]);
                        Hy[k][j] = (dNdy[j] * dNdy[k]);
                        H[k][j] = (Hx[k][j] + Hy[k][j]);
                    }
                }
                System.out.println("Matrix Hx:");
                showMatrix(Hx);
                System.out.println("Matrix Hy:");
                showMatrix(Hy);
                System.out.println("Matrix Hx+Hy:");
                showMatrix(H);

                System.out.println("Matrix H | k*(Hx+Hy)*detJ");
                for (int k = 0; k < 4; k++) {
                    for (int j = 0; j < 4; j++) {
                        Hlocal[k][j] = (wspK * H[k][j] * det);
                    }
                }
                showMatrix(Hlocal);

                for (int k = 0; k < 4; k++) {
                    for (int j = 0; j < 4; j++) {
                        Hsum[k][j] += Hlocal[k][j];
                        //System.out.println("Hsum[" + k + "][" + j + "] = " + Hsum[k][j]);
                    }
                }

            }

            int[] BC = elements[i-1].getArea();
            if(BC[0] == 1) {
                shapeFunction(-1, (1/sqrt(3)), alfa, delta);
                shapeFunction(-1, -(1/sqrt(3)), alfa, delta);
            }
            if(BC[1] == 1) {
                shapeFunction(-(1/sqrt(3)), -1, alfa, delta);
                shapeFunction((1/sqrt(3)), -1, alfa, delta);
            }
            if(BC[2] == 1) {
                shapeFunction(1, -(1/sqrt(3)), alfa, delta);
                shapeFunction(1, (1/sqrt(3)), alfa, delta);
            }
            if(BC[3] == 1) {
                shapeFunction(-(1/sqrt(3)), 1, alfa, delta);
                shapeFunction((1/sqrt(3)), 1, alfa, delta);
            }
            System.out.println("Integral for each boundary condition | alfa*{N}*{N}^T*(delta/2):");
            showMatrix(Nsum);

            System.out.println("Matrix H for the " + i + " element (sum from all integration points):");
            showMatrix(Hsum);

            for (int k = 0; k < 4; k++) {
                for (int j = 0; j < 4; j++) {
                    Hsum[k][j] = (Hsum[k][j] + Nsum[k][j]);
                }
            }

            System.out.println("Local matrix H (final for the " + i + " element):");
            showMatrix(Hsum);

            int [] elementID = elements[i-1].getID();

            for (int k = 0; k < 4; k++) {
                for (int j = 0; j < 4; j++) {
                    int id_k = nodes[elementID[k] - 1].getI();
                    int id_j = nodes[elementID[j] - 1].getI();
                    Hglobal[id_k-1][id_j-1] += Hsum[k][j];
                }
            }

            for (int k = 0; k < 4; k++) {
                for (int j = 0; j < 4; j++) {
                    Hsum[k][j] = 0;
                    Nsum[k][j] = 0;
                }
            }

        }
        System.out.println("Global matrix [H] :");
        for (int k = 0; k < Hglobal.length; k++) {
            for (int j = 0; j < Hglobal.length; j++) {
                System.out.printf("%.4f   ", Hglobal[k][j]);
            }
            System.out.println();
        }
        System.out.println();
    }

    public void matrixC(Grid grid){
        double[] ksi = {-(1/sqrt(3)), (1/sqrt(3)), (1/sqrt(3)), -(1/sqrt(3))};
        double[] eta = {-(1/sqrt(3)), -(1/sqrt(3)), (1/sqrt(3)), (1/sqrt(3))};
        double[][] C = new double[4][4];
        double[][] Csum = new double[4][4];

        int cieplo_wlasciwe = grid.getGlobalData().getCieplo_wlasciwe();
        int gestosc = grid.getGlobalData().getGestosc();

        for(int i = 1; i <= elements.length; i++) {
            double[] shapeF = new double[4];
            for (int p = 0; p < 4; p++) {
                shapeF[0] = 0.25 * (1 - ksi[p]) * (1 - eta[p]);
                shapeF[1] = 0.25 * (1 + ksi[p]) * (1 - eta[p]);
                shapeF[2] = 0.25 * (1 + ksi[p]) * (1 + eta[p]);
                shapeF[3] = 0.25 * (1 - ksi[p]) * (1 + eta[p]);


                for (int k = 0; k < 4; k++) {
                    for (int j = 0; j < 4; j++) {
                        C[k][j] += shapeF[k] * shapeF[j] * det * cieplo_wlasciwe * gestosc;
                    }
                }

                for (int k = 0; k < 4; k++) {
                    for (int j = 0; j < 4; j++) {
                        Csum[k][j] += C[k][j];
                    }
                }
                resetMatrix(C);
            }

            int [] elementID = elements[i-1].getID();

            for (int k = 0; k < 4; k++) {
                for (int j = 0; j < 4; j++) {
                    int id_k = nodes[elementID[k] - 1].getI();
                    int id_j = nodes[elementID[j] - 1].getI();
                    Cglobal[id_k-1][id_j-1] += Csum[k][j];
                }
            }
            resetMatrix(Csum);
        }
        System.out.println("Global matrix [C] :");
        for (int k = 0; k < Cglobal.length; k++) {
            for (int j = 0; j < Cglobal.length; j++) {
                System.out.printf("%.4f   ", Cglobal[k][j]);
            }
            System.out.println();
        }
        System.out.println();

    }

    public void vectorP(Grid grid){
        double alfa = grid.getGlobalData().getAlfa();
        double B = grid.getGlobalData().getB();
        int nB = grid.getGlobalData().get_nB() - 1; //liczba wezlow na szerokosc -1
        double delta = B/nB; //długosc boku
        this.dt = grid.getGlobalData().getDt();

        int tempOtoczenia = grid.getGlobalData().getTempOtoczenia();

        for(int i = 1; i <= elements.length; i++) {
            System.out.println("ELEMENT " + i);
            int[] BC = elements[i - 1].getArea();
            if (BC[0] == 1) {
                System.out.println("Sciana 1");
                shapeFunction2(-1, 1 / sqrt(3), 0);
                shapeFunction2(-1, -1 / sqrt(3), 0);
            }
            if (BC[1] == 1) {
                System.out.println("Sciana 2");
                shapeFunction2(-(1 / sqrt(3)), -1, 1);
                shapeFunction2((1 / sqrt(3)), -1, 1);
            }
            if (BC[2] == 1) {
                System.out.println("Sciana 3");
                shapeFunction2(1, -(1 / sqrt(3)), 2);
                shapeFunction2(1, (1 / sqrt(3)), 2);
            }
            if (BC[3] == 1) {
                System.out.println("Sciana 4");
                shapeFunction2(-(1 / sqrt(3)), 1, 3);
                shapeFunction2((1 / sqrt(3)), 1, 3);
            }
            System.out.println("warBrzeg[][]:");
            for(int k = 0; k < 4; k++) {
                for (int j = 0; j < 4; j++) {
                    System.out.print(warBrzeg[k][j] + "   ");
                }
                System.out.println();
            }
            System.out.println();

            for(int k = 0; k < 4; k++){
                for(int j = 0; j < 4; j++){
                    P[k] += (alfa * tempOtoczenia * (delta/2) * warBrzeg[k][j]);
                }
            }

            System.out.println("Wektor lokalny P:");
            for(int j = 0; j < 4; j++){
                System.out.println(P[j]);
            }
            System.out.println();

            int [] elementID = elements[i-1].getID();
            for (int j = 0; j < 4; j++) {
                    int id_j = nodes[elementID[j] - 1].getI();
                    Pglobal[id_j-1] += P[j];
            }

            resetMatrix(warBrzeg);
            for(int j = 0; j < 4; j++){
                P[j] = 0;
            }
        }
        System.out.println("Global vector {P} :");
        for (int j = 0; j < Pglobal.length; j++) {
            System.out.print(Pglobal[j] + "   ");
        }
        System.out.println("\n");
    }

    public void finalMatrixH(){
        int dt = 50;

        for(int i = 0; i < ND; i++){
            for(int j = 0; j < ND; j++){
                Hfinal[i][j] = Hglobal[i][j] + Cglobal[i][j]/dt;
            }
        }
        System.out.println("Final matrix [H] = [H]+[C]/dt :");
        for(int i = 0; i < ND; i++){
            for(int j = 0; j < ND; j++){
                System.out.print(Hfinal[i][j] + "   ");
            }
            System.out.println();
        }
    }

    public void finalVectorP(){
        int time = 1;
        double[] Pfinal;
        T0 = new double[ND];
        for (int k = 0; k < ND; k++) {
                T0[k] = 25.;
        }

        //System.out.println("Final vector {P} = {P}+([C]/dt)*{T0} - iteration 0:");
        Pfinal = calculateP(T0);

        for(int k = 0; k < 30; k ++){
            double[] T1 = Gauss.gaussEliminationMethod(Hfinal, Pfinal, ND);

            double minT = T1[0];
            double maxT = T1[0];
            for(int i = 0; i < ND; i++) {
                if(T1[i] < minT){
                    minT = T1[i];
                }
                if(T1[i] > maxT){
                    maxT = T1[i];
                }
            }

            System.out.println("Time: " + time + " min.");
            System.out.println("minT: " + minT);
            System.out.println("maxT: " + maxT);

            //System.out.print(minT + "\t");
            //System.out.print(maxT);

            Pfinal = calculateP(T1);

            time += 1;
        }

    }

    private double[] calculateP(double[] vectorT){
        double[] Pfinal = new double[ND];

        for(int i = 0; i < ND; i++) {
            for (int j = 0; j < ND; j++) {
                Pfinal[i] += (Cglobal[i][j]/dt) * vectorT[j];
            }
        }
        for(int i = 0; i < ND; i++) {
            Pfinal[i] += Pglobal[i];
        }
        System.out.println();

        return Pfinal;
    }

    private void shapeFunction2(double ksi, double eta, int i) {
        double[] shapeF = new double[4];
        shapeF[0] = 0.25 * (1 - ksi) * (1 - eta);
        shapeF[1] = 0.25 * (1 + ksi) * (1 - eta);
        shapeF[2] = 0.25 * (1 + ksi) * (1 + eta);
        shapeF[3] = 0.25 * (1 - ksi) * (1 + eta);

        for(int j = 0; j < 4; j++){
            warBrzeg[i][j] += shapeF[j];
        }
    }

    private void shapeFunction(double ksi, double eta, double alfa, double delta) {
        double[] shapeF = new double[4];
        shapeF[0] = 0.25 * (1 - ksi) * (1 - eta);
        shapeF[1] = 0.25 * (1 + ksi) * (1 - eta);
        shapeF[2] = 0.25 * (1 + ksi) * (1 + eta);
        shapeF[3] = 0.25 * (1 - ksi) * (1 + eta);

        for (int k = 0; k < 4; k++) {
            for (int j = 0; j < 4; j++) {
                N[k][j] += (shapeF[k] * shapeF[j] * (delta / 2) * alfa);
            }
        }
        for (int k = 0; k < 4; k++) {
            for (int j = 0; j < 4; j++) {
                Nsum[k][j] += N[k][j];
            }
        }
        resetMatrix(N);
    }

    private void resetMatrix(double [][] tab){
        for (int k = 0; k < 4; k++) {
            for (int j = 0; j < 4; j++) {
                tab[k][j] = 0.;
            }
        }
    }

    private void showMatrix(double [][] tab){
        for(int k = 0; k < 4; k++) {
            for (int j = 0; j < 4; j++) {
                System.out.printf("%.6f   ", tab[k][j]);
            }
            System.out.println();
        }
        System.out.println();
    }

    private void showJacobian(int i, int p, double det, double[][][][] inverseMatrix){
        System.out.println("the integration point: " + (p+1));
        System.out.printf("%.6f   ", J[i-1][p][0][0]);
        System.out.printf("%.6f \n", J[i-1][p][0][1]);
        System.out.printf("%.6f   ", J[i-1][p][1][0]);
        System.out.printf("%.6f \n", J[i-1][p][1][1]);

        System.out.println("det of Jacobian: " + det);

        System.out.println("Inverse matrix:");
        System.out.printf("%.6f   ", inverseMatrix[i-1][p][0][0]);
        System.out.printf("%.6f \n", inverseMatrix[i-1][p][0][1]);
        System.out.printf("%.6f   ", inverseMatrix[i-1][p][1][0]);
        System.out.printf("%.6f \n", inverseMatrix[i-1][p][1][1]);
    }
}
