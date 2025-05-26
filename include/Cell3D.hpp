#include <vector>
#include <array>
#include <array>
#include <glm/vec3.hpp>
#include <glm/vec2.hpp>
#include "ECM.hpp"

#ifndef __CELL3D__
#define __CELL3D__
namespace DPM3D {
    class Cell
    {
        public:
            int NV;
            double calA0;
            double v0;
            double r0;
            double a0;
            double s0;
            double l0;
            double U;
            double Kv;
            double Ka;
            double Kl;
            double Kb;
            double Ks;
            int ntriangles;
            std::vector<std::vector<int>> midpointCache;
            std::vector<glm::dvec3> Positions;
            std::vector<glm::dvec3> Velocities;
            std::vector<glm::dvec3> Forces;
            std::vector<std::vector<int>> FaceIndices;
            std::vector<bool> isJunction;
            std::vector<bool> isFocalAdh;

            Cell(double x1,double y1, double z1,
                double calA0,int f,double r0,
                double Kv, double Ka,double Kb);

            Cell(std::vector<double>X,
                 std::vector<double>Y,
                 std::vector<double>Z,
                 std::vector<std::vector<int>> Triangles);

            void ResetForces();
            void AddFaceIndex(int a, int b, int c);
            glm::dvec3 GetMiddlePoint(int i, int j);
            int AddMiddlePoint(int p1, int p2);
            void AddVertex(glm::dvec3 vec);
            void ShapeForceUpdate();
            void VolumeForceUpdate();
            void AreaForceUpdate();
            void BendingForceUpdate();
            void FindJunctions();
            std::array<std::vector<double>,2> getJunctionPoints();
            void FindFocalAdhesion();
            void SurfaceGradient(double x, double mindist, double G);
            void StickToSurface(double x, double mindist);
            void StickToSurface(DPM3D::ECM ECM, double mindist);
            void StickToSurfaceCatch(double z, double mindist);
            void StickToSurfaceCatch(DPM3D::ECM ECM, double mindist);
            void StickToSurfaceSlip(double z, double mindist);
            void StickToSurfaceSlip(DPM3D::ECM ECM, double mindist);
            void SetupSurface(int n);
            void RepelSurface(double z);
            void SurfaceStrech(double scale);
            void ExtendVertex(int vi, double force);
            void EulerUpdate(int steps, double dt);
            void EulerUpdate(double dt);
            void FIREMinimization(double alpha, double dt, int itmax, double Ftol);
            double GetVolume();
            double GetArea(int triidx);
            double GetSurfaceArea();
            double GetCalA0();
            bool pointInside(glm::dvec3);
            glm::dvec3 GetCOM();
            static std::vector<glm::dvec3> GetShapeForces(DPM3D::Cell Cell);
            static std::vector<glm::dvec3> GetVolumeForces(DPM3D::Cell Cell);
            static std::vector<glm::dvec3> GetAreaForces(DPM3D::Cell Cell);
            static std::vector<glm::dvec3> GetBendingForces(DPM3D::Cell Cell);
            std::array<std::vector<double>,3> GetPositions();
            std::array<std::vector<double>,3> GetForces();
    };

    glm::dvec3 GetMiddlePoint(glm::dvec3 a, glm::dvec3 b);
    bool polarCompare(glm::dvec2& a, glm::dvec2& b, glm::dvec2& p0);
    std::vector<int> newFace(int a, int b, int c);
    double urandom();
    double norm(glm::dvec3 vec);
    double norm2(glm::dvec3 vec);
    double distance(glm::dvec3 a,glm::dvec3 b);
    double distanceSq(glm::dvec3 a,glm::dvec3 b);
    int orientation(glm::dvec2 a, glm::dvec2 b, glm::dvec2 c);

}

#endif
