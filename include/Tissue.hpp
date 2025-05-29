#ifndef __TISSUE__
#define __TISSUE__

#include <vector>
#include "Cell3D.hpp"
#include "ECM.hpp"
#include <string>
#include <array>

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
        std::string attactionMethod;
        std::vector<DPM3D::Cell> Cells;

        Tissue(std::vector<DPM3D::Cell> Cells, double phi0);

        void CellInteractingUpdate(int ci);
        void UpdateJunctions();
        void InteractingUpdate();
        void UpdateShapeForces();
        void EulerUpdate(int nsteps, double dt);
        void EulerUpdate(double dt);
        std::vector<bool> FindOverlaps(int ci, int cj);
        std::vector<double> findWindingNumber(int ci, int cj);
        void MonolayerDisperse();
        void TissueDisperse();
        void StickToSurface(double z, double mindist);
        void CellStickToSurface(int ci, double z, double mindist);
        void InteractECM(DPM3D::ECM ecm);
        void JunctionSlipForceUpdate(int ci);
        void JunctionCatchForceUpdate(int ci);
        void GeneralAttraction(int ci);
        std::vector<std::vector<double>> GetVesselPosition(int ci);
        void applyShearStress(std::array<double,3> fluidVelocity, double viscosity, double velocityGradient, double stiffness);
    };
}

#endif
