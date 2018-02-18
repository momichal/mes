import java.util.Arrays;

public class Element {
    private int[] ID = new int[4];
    private int[] area = new int[4];    // powierzchnia BC
    private int i;

    public Element(int[] ID, int[] area, int i) {
        this.ID = ID;
        this.area = area;
        this.i = ++i;
    }

    public int[] getID() {
        return ID;
    }

    public int[] getArea() {
        return area;
    }

    @Override
    public String toString() {
        return "El. " + String.format("%d", i) +
                " {" +
                "ID = " + Arrays.toString(ID) + "   " +
                "BC ścian = " + Arrays.toString(area) +
                '}';
    }
}
