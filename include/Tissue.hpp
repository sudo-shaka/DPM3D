#include <vector>
#include <Cell3D.hpp>

namespace DPM3D{
    class Tissue{
        public:
        double phi0;
        int NCELLS;
        int VertDOF;
        double L;
        double Kc;
        double U;
        std::vector<DPM3D::Cell> Cells;

        Tissue(std::vector<DPM3D::Cell> Cells, double phi0);

        void InteractingForceUpdate();
    };
}