import java.io.File;
import java.io.IOException;
import java.util.Scanner;

public class GlobalData {
    private  double H, B;   // Wysokość, szerokość siatki
    private int nH, nB;     // Liczba węzłów na wysokość, szerokość
    private  int nh, ne;    // Liczba węzłów, elementów
    private double t;
    private double k;       // Przewodność
    private double alfa;
    private int dt;
    private int cieplo_wlasciwe;
    private int gestosc;
    private int tempOtoczenia;

    public GlobalData(String filePath) throws IOException {
        File file = new File(filePath);
        try (Scanner scanner = new Scanner(file)) {
            H = scanner.nextDouble(); scanner.nextLine();
            B = scanner.nextDouble(); scanner.nextLine();
            nH = scanner.nextInt(); scanner.nextLine();
            nB = scanner.nextInt(); scanner.nextLine();
            t = scanner.nextDouble(); scanner.nextLine();
            k = scanner.nextDouble(); scanner.nextLine();
            alfa = scanner.nextDouble(); scanner.nextLine();
            dt = scanner.nextInt(); scanner.nextLine();
            cieplo_wlasciwe = scanner.nextInt(); scanner.nextLine();
            gestosc = scanner.nextInt(); scanner.nextLine();
            tempOtoczenia = scanner.nextInt();
        }
        nh = nH * nB;
        ne = (nH - 1) * (nB - 1);
    }

    public double getTemp() {
        return t;
    }

    public  double getH() {
        return H;
    }

    public  double getB() {
        return B;
    }

    public  int get_nH() {
        return nH;
    }

    public  int get_nB() {
        return nB;
    }

    public  int get_nh() {
        return nh;
    }

    public  int get_ne() {
        return ne;
    }

    public double getK() {
        return k;
    }

    public double getAlfa() {
        return alfa;
    }

    public int getDt() {
        return dt;
    }

    public int getCieplo_wlasciwe() {
        return cieplo_wlasciwe;
    }

    public int getGestosc() {
        return gestosc;
    }

    public int getTempOtoczenia() {
        return tempOtoczenia;
    }
}
