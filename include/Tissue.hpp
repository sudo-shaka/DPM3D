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
        double phi0;  //initial packing fraction
        int NCELLS;   //number of cells
        double L;     //boxlength for periodic boundary conditions (PBC)
        double Kre;   //Repulsion stiffness
        double Kat;   //Attractive stiffness
        bool PBC = true; //if PBC will be applied
        double f0 = 1.0; //charateristic force scale for bond behavoir
        std::string attactionMethod; //Method of attraction between verts
        std::vector<DPM3D::Cell> Cells; //List of cell objects that make up the tissue

        Tissue(std::vector<DPM3D::Cell> Cells, double phi0); //constructor

        void CellInteractingUpdate(int ci); //Force update for intercellular interactions for a given cell (ci)
        void UpdateJunctions();   //Find the verts which make the perimeter looking down in Z direciton
        void InteractingUpdate(); //Force update for intercelluar interactions for all cells (additive)
        void UpdateShapeForces(); //Get shape force update for all cells (additive)
        void EulerUpdate(int nsteps, double dt); //Updates positions based on force, resets forces
        void EulerUpdate(double dt); // same as above but only once
        std::vector<bool> FindOverlaps(int ci, int cj); //get a list of bools indicating a vertex in all verts (in order) overlaps with another cell. Shorthand : OverlappingVerts = AllVerts[overlaps]
        std::vector<double> findWindingNumber(int ci, int cj); //This finds the winding number for all verts in ci for overlapping with cells cj. winding number affects the force of Repulsion as well as determines of vert is overlapping with another cell
        void MonolayerDisperse(); //randomly distributes all cells (minimizing overlaps) in the X/Y plane
        void TissueDisperse(); //distributes all cells over X/Y/Z plane
        void StickToSurface(double z, double mindist); //z is the distance of surface, mindist is the minimum distance for attraction to occour
        void CellStickToSurface(int ci, double z, double mindist); // same as above for a specific cell in tissue
        void InteractECM(DPM3D::ECM ecm); //ECM is a list a point for a cells to adhere to
        void JunctionSlipForceUpdate(int ci); //Attraction method for cell (ci)
        void JunctionCatchForceUpdate(int ci); ////Attraction method for cell (ci)
        void GeneralAttraction(int ci);//Attraction method for cell (ci)
        std::vector<std::vector<double>> GetVesselPosition(int ci); //Takes cells attached to z=0 surface with PBC to make a cylindircal surface of cells (plotting plane as a vessel)
        void applyShearStress(std::array<double,3> fluidVelocity, double viscosity, double velocityGradient, double stiffness); //Apply a shear force over the faces of the cells going from the x direction
    };
}

#endif
