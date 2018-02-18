import static java.lang.StrictMath.sqrt;

public class UniversalElement {
    private double[][] shapeFunc;
    private double[] ksi = {-(1/sqrt(3)), (1/sqrt(3)), (1/sqrt(3)), -(1/sqrt(3))};
    private double[] eta = {-(1/sqrt(3)), -(1/sqrt(3)), (1/sqrt(3)), (1/sqrt(3))};
    private double[][] derivativeKsi;
    private double[][] derivativeEta;

    public UniversalElement(){
        shapeFunc = new double[4][4];
        derivativeKsi = new double[4][4];
        derivativeEta = new double[4][4];
    }

    public void shapeFunc() {
        for(int i=0; i<4; i++) {
            shapeFunc[i][0] = 0.25 * (1 - ksi[i]) * (1 - eta[i]);
            shapeFunc[i][1] = 0.25 * (1 + ksi[i]) * (1 - eta[i]);
            shapeFunc[i][2] = 0.25 * (1 + ksi[i]) * (1 + eta[i]);
            shapeFunc[i][3] = 0.25 * (1 - ksi[i]) * (1 + eta[i]);
        }
    }

    public void derivativeKsi() {
        for(int i=0; i<4; i++) {
            derivativeKsi[i][0] = -0.25 * (1 - eta[i]);
            derivativeKsi[i][1] = 0.25 * (1 - eta[i]);
            derivativeKsi[i][2] = 0.25 * (1 + eta[i]);
            derivativeKsi[i][3] = -0.25 * (1 + eta[i]);
        }
    }

    public void derivativeEta(){
        for(int i=0; i<4; i++) {
            derivativeEta[i][0] = -0.25 * (1 - ksi[i]);
            derivativeEta[i][1] = -0.25 * (1 + ksi[i]);
            derivativeEta[i][2] = 0.25 * (1 + ksi[i]);
            derivativeEta[i][3] = 0.25 * (1 - ksi[i]);
        }
    }

    public void showShapeFunc() {
        System.out.println("Funkcje kształtu:");
        for (int i = 0; i < 4; i++) {
            System.out.print("Pkt " + (i+1) + ": ");
            for (int j = 0; j < 4; j++)
                System.out.printf("%.6f   ", shapeFunc[i][j]);
        System.out.println();
        }
    }

    public void showDerivativeKsi(){
        System.out.println("\nPochodne po ksi:");
        System.out.println("\t\t N1 \t\t N2 \t\t N3 \t\t N4");
        for (int i = 0; i < 4; i++) {
            System.out.print("Pkt " + (i+1) + ": ");
            for (int j = 0; j < 4; j++)
                System.out.printf("%.6f   ", derivativeKsi[i][j]);
            System.out.println();
        }
    }

    public void showDerivativeEta(){
        System.out.println("\nPochodne po eta:");
        System.out.println("\t\t N1 \t\t N2 \t\t N3 \t\t N4");
        for (int i = 0; i < 4; i++) {
            System.out.print("Pkt " + (i+1) + ": ");
            for (int j = 0; j < 4; j++)
                System.out.printf("%.6f   ", derivativeEta[i][j]);
            System.out.println();
        }
    }

    public double[][] getDerivativeKsi() {
        return derivativeKsi;
    }

    public double[][] getDerivativeEta() {
        return derivativeEta;
    }
}
