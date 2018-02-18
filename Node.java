public class Node {
    private double x, y;
    private double t;
    private boolean boundaryCondition; //warunek brzegowy
    private int i;

    public Node(double x, double y, double t, boolean boundaryCondition, int i) {
        this.x = x;
        this.y = y;
        this.t = t;
        this.boundaryCondition = boundaryCondition;
        this.i = ++i;
    }

    public double getX() {
        return x;
    }

    public double getY() {
        return y;
    }

    public int getI() {
        return i;
    }

    @Override
    public String toString() {
        return  "Node " + String.format("%d", i) +
                " {" +
                "x=" + String.format("%.4f", x) +
                ", y=" + String.format("%.4f", y) +
                ", t=" + String.format("%.4f", t) +
                ", BC=" + boundaryCondition +
                '}';
    }
}
