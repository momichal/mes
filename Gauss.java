public class Gauss {
    private static final double e = Math.pow(10, -12);

    public static double[] gaussEliminationMethod(double[][] A, double[] B, int size) {
        double[] result = new double[size];
        double[][] U = new double[size][size + 1];


        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                U[j][i] = A[j][i];
            }
        }

        for (int i = 0; i < size; ++i) {
            U[i][size] = B[i];
        }

        for (int i = 0; i < size - 1; ++i) {
            for (int j = i + 1; j < size; ++j) {
                if (Math.abs(U[i][i]) < e) {
                    break;
                }

                double x = -U[j][i] / U[i][i];
                for (int k = 0; k < size + 1; ++k) {
                    U[j][k] += x * U[i][k];
                }
            }
        }

        for (int i = size - 1; i >= 0; --i) {
            double x = U[i][size];
            for (int j = size - 1; j >= 0; --j) {
                x -= U[i][j] * result[j];
            }
            if(Math.abs(U[i][i]) < e) {
                break;
            }
            result[i] = x / U[i][i];
        }

        return result;
    }

}
