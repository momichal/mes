import java.io.IOException;

public class Main
{
    public static void main(String[] args)
    {
        GlobalData globalData = null;
        try {
            globalData = new GlobalData("src/data.txt");
        } catch (IOException e) {
            e.printStackTrace();
        }

        /* Creating a grid */
        Grid grid = new Grid(globalData);
        grid.initializeNodes();
        grid.initializeElements();

        grid.viewNodes();
        System.out.println();
        grid.viewElements();
        System.out.println();

        /* Shape function */
        UniversalElement universalElement = new UniversalElement();
        universalElement.shapeFunc();
        universalElement.showShapeFunc();

        /* Derivatives after ksi, eta */
        universalElement.derivativeKsi();
        universalElement.derivativeEta();
        universalElement.showDerivativeKsi();
        universalElement.showDerivativeEta();
        System.out.println();

        /* Jacobian */
        Jacobian jacobian = new Jacobian(grid);
        jacobian.jacobian(universalElement, grid);

        /* Matrix H */
        jacobian.transformationJacobian(grid);

        /* Matrix C */
        jacobian.matrixC(grid);

        /* Vector P */
        jacobian.vectorP(grid);

        jacobian.finalMatrixH();
        jacobian.finalVectorP();
    }
}
