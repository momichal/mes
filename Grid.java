public class Grid {
    private GlobalData globalData;
    private Node[] ND;
    private Element[] EL;

    public Grid(GlobalData globalData) {
        this.globalData = globalData;
        ND = new Node[globalData.get_nh()];
        EL = new Element[globalData.get_ne()];
    }

    public void initializeNodes() {
        double x = 0, y = 0;
        double dx = globalData.getB() / (globalData.get_nB() - 1);
        double dy = globalData.getH() / (globalData.get_nH() - 1);

        for (int i = 0; i < ND.length; ) {
            ND[i] = new Node(x, y, globalData.getTemp(), statusBoundaryCondition(i), i);
            y += dy;
            ++i;
            if (i % globalData.get_nH() == 0) {
                x += dx;
                y = 0;
            }
        }
    }

    private boolean statusBoundaryCondition(int i) {
        boolean BC = false;

        if (i < globalData.get_nH() || i >= globalData.get_nh() - globalData.get_nH()) {
            BC = true;
        }
        if (i % globalData.get_nH() == 0 || i % globalData.get_nH() == globalData.get_nH() - 1) {
            BC = true;
        }
        return BC;
    }

    public void initializeElements() {
        int idB = 0, idH = 0;
        for (int i = 0; i < EL.length; i++) {
            if (i % (globalData.get_nH() - 1) == 0) {
                idH = 0;
                ++idB;
            }
            ++idH;
            EL[i] = new Element(elementID(idH, idB), elementBC(idH, idB), i);
        }
    }

    private int[] elementID(int idH, int idB) {

        int nH = globalData.get_nH();
        int down_left = nH * (idB - 1) + idH;
        int down_right = down_left + nH;
        int up_right = down_right + 1;
        int up_left = down_left + 1;
        return new int[]{down_left, down_right, up_right, up_left};
    }

    private int[] elementBC(int idH, int idB) {
        int[] BC = new int[]{0, 0, 0, 0}; //ściany
        if (idB == 1) {
            BC[0] = 1;
        }
        if (idH == 1) {
            BC[1] = 1;
        }
        if (idB == (globalData.get_nB() - 1)) {
            BC[2] = 1;
        }
        if (idH == (globalData.get_nH() - 1)) {
            BC[3] = 1;
        }
        return BC;
    }

    public void viewNodes() {
        for (int i = 0; i < ND.length; i++)
            System.out.println(ND[i]);
    }

    public void viewElements() {
        for (int i = 0; i < EL.length; i++)
            System.out.println(EL[i]);
    }

    public GlobalData getGlobalData() {
        return globalData;
    }

    public Node[] getND() {
        return ND;
    }

    public Element[] getEL() {
        return EL;
    }
}
