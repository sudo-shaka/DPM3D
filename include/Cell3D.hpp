#include <vector>
#include <glm/vec3.hpp>
namespace DPM3D {
    class Cell
    {
        public:
            int NV;
            double calA0;
            double v0;
            double a0;
            double s0;
            double U;
            double Kv;
            double Ka;
            double Kl;
            double Kb;
            int ntriangles;
            std::vector<std::vector<int>> midpointCache;
            std::vector<glm::dvec3> Positions;
            std::vector<glm::dvec3> Velocities;
            std::vector<glm::dvec3> Forces;
            std::vector<std::vector<int>> FaceIndices;
            std::vector<glm::dvec3> surfacepositions;
            int  nsurfacep;


            Cell(double x1,double y1, double z1,
                double calA0,int f,
                double Kv, double Ka,double Kb);

            void ResetForces();
            void AddFaceIndex(int a, int b, int c);
            glm::dvec3 GetMiddlePoint(int i, int j);
            int AddMiddlePoint(int p1, int p2);
            void AddVertex(glm::dvec3 vec);
            void ShapeForceUpdate();
            void VolumeForceUpdate();
            void AreaForceUpdate();
            void BendingForceUpdate();
            void StickToSurface(double x, double mindist);
            void SetupSurface(double z);
            void SurfaceStrech(double scale);
            void ExtendVertex(int vi, double force);
            void Crawling();
            void EulerUpdate(int steps, double dt);
            void EulerUpdate(double dt);
            void FIREMinimization(double alpha, double dt, int itmax, double Ftol);
            double GetVolume();
            double GetArea(int triidx);
            double GetSurfaceArea();
            double GetCalA0();
            glm::dvec3 GetCOM();

    };
    glm::dvec3 GetMiddlePoint(glm::dvec3 a, glm::dvec3 b);
    std::vector<int> newFace(int a, int b, int c);
    double urandom();
    double norm(glm::dvec3 vec);
    double norm2(glm::dvec3 vec);
    double distance(glm::dvec3 a,glm::dvec3 b);
    double distanceSq(glm::dvec3 a,glm::dvec3 b);

}
