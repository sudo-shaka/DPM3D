#include <vector>
#include <Cell3D.hpp>
#include <ECM.hpp>

namespace DPM3D{
    class Tissue{
        public:
        double phi0;
        int NCELLS;
        int VertDOF;
        double L;
        double Kre;
        double Kat;
        double U;
        bool PBC;
        std::vector<DPM3D::Cell> Cells;

        Tissue(std::vector<DPM3D::Cell> Cells, double phi0);

        void CellInteractingUpdate(int ci);
        void InteractingUpdate();
        void UpdateShapeForces();
        void EulerUpdate(int nsteps, double dt);
        void EulerUpdate(double dt);
        std::vector<bool> FindOverlaps(int ci, int cj);
        void MonolayerDisperse();
        void TissueDisperse();
        void StickToSurface(double z, double mindist);
        void CellStickToSurface(int ci, double z, double mindist);
        void InteractECM(DPM3D::ECM ecm);
    };
}
