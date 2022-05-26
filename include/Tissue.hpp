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
        std::vector<bool> overlaps;

        Tissue(std::vector<DPM3D::Cell> Cells, double phi0);

        void RetractingForceUpdate();
        void CellRetractingUpdate(int ci);
        void UpdateShapeForces();
        void EulerUpdate(int nsteps, double dt);
        void EulerUpdate(double dt);
        void FindOverlaps(int ci, int cj);
        void MonolayerDisperse();
        void TissueDisperse();
    };
}