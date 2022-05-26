#include <cmath>
#include <Tissue.hpp>
#include <vector>
#include <iostream>
#include <glm/glm.hpp>
#include <glm/vec3.hpp>
#include <glm/gtx/norm.hpp>
#include <thread>

namespace DPM3D{
    Tissue::Tissue(std::vector<DPM3D::Cell> _cells, double _phi0){
        Cells = _cells;
        phi0 = _phi0;
        NCELLS = Cells.size();
        double volume = 0.0;
        for(int i=0; i<NCELLS;i++){
            volume+=(4/3)*M_PI*pow(Cells[i].r0,3);
        }
        L = cbrt(volume)/phi0;
    }

    void Tissue::MonolayerDisperse(){
        std::vector<double> X,Y,Fx,Fy;
        X.resize(NCELLS); Y.resize(NCELLS);
        Fx.resize(NCELLS); Fy.resize(NCELLS);
        double ri,rj,yi,yj,xi,xj,dx,dy,dist;
        double ux,uy,ftmp,fx,fy;
        int  i,j,count;
        for(i=0;i<NCELLS;i++){
            X[i] = drand48() * L;
            Y[i] = drand48() * L;
        }
        double oldU = 100, dU = 100;
        count = 0;
        while(dU > 1e-6){
            U = 0;
            for(i=0;i<NCELLS;i++){
                Fx[i] = 0.0;
                Fy[i] = 0.0;
            }

            for(i=0;i<NCELLS;i++){
                xi = X[i];
                yi = Y[i];
                ri = Cells[i].r0;
                for(j=0;j<NCELLS;j++){
                    if(j != i){
                        xj = X[j];
                        yj = Y[j];
                        rj = Cells[j].r0;
                        dx = xj-xi;
                        dx -= L*round(dx/L);
                        dy = yj-yi;
                        dy -= L*round(dy/L);
                        dist = sqrt(dx*dx + dy*dy);
                        if(dist < 0.0)
                            dist *= -1;
                        if(dist <= (ri+rj)){
                            ux = dx/dist;
                            uy = dy/dist;
                            ftmp = (1.0-dist/(ri+rj))/(ri+rj);
                            fx = ftmp*ux;
                            fy = ftmp*uy;
                            Fx[i] -= fx;
                            Fy[i] -= fy;
                            Fy[j] += fy;
                            Fx[j] += fx;
                            U += 0.5*(1-(dist/(ri+rj))*(1-dist/(ri+rj)));
                        }
                    }
                }
            }
            for(int i=0; i<NCELLS;i++){
                X[i] += 0.01*Fx[i];
                Y[i] += 0.01*Fy[i];
            }
            dU = U-oldU;
            if(dU < 0.0)
                dU *= -1;
            oldU = U;
            count++;
            if(count > 1e5){
                std::cerr << "Warning: Max timesteps for dispersion reached"  << std::endl;
                break;
            }
        }
        glm::dvec3 com;
        for(i=0; i<NCELLS; i++){
            com = Cells[i].GetCOM();
            for(j=0;j<Cells[i].NV;j++){
                Cells[i].Positions[j].x -= com.x;
                Cells[i].Positions[j].y -= com.y;
                Cells[i].Positions[j].x += X[i];
                Cells[i].Positions[j].y += Y[i];
            }
        }
    }

    void Tissue::TissueDisperse(){
        std::vector<glm::dvec3> centers;
        std::vector<glm::dvec3> forces;
        glm::dvec3 rij;
        centers.resize(NCELLS);
        forces.resize(NCELLS);
        int  i,j,count=0;
        double ftmp;
        for(i=0;i<NCELLS;i++){
            centers[i].x = drand48() * L;
            centers[i].y = drand48() * L;
            centers[i].z = drand48() * L;
        }
        double oldU = 100, dU = 100, U, dist;
        while(dU > 1e-6){
            U = 0;
            for(i=0;i<NCELLS;i++){
                forces[i] = {0,0,0};
            }
            for(i=0;i<NCELLS;i++){
                for(j=0;j<NCELLS;j++){
                    if(j!=i){
                        rij = centers[j] - centers[i];
                        rij -= L*round(rij/L);
                        dist = sqrt(rij.x*rij.x + rij.y*rij.y + rij.z*rij.z);
                        if(dist < 0.0){
                            dist *=-1;
                        }
                        if(dist < (Cells[i].r0 + Cells[j].r0)){
                            ftmp = (1-dist/(Cells[i].r0+Cells[j].r0)/(Cells[i].r0+Cells[j].r0));
                            forces[i] -= ftmp*glm::normalize(rij);
                            forces[j] += ftmp*glm::normalize(rij);
                            U += 0.5*(1-(dist/(Cells[i].r0+Cells[j].r0))*(1-dist/(Cells[i].r0+Cells[j].r0)));
                        }
                    }
                }
            }
            for(i=0;i<NCELLS;i++){
                centers[i] += 0.01*forces[i];
            }
            dU = U-oldU;
            if(dU < 0.0){
                dU *= -1;
            }
            oldU = U;
            count++;
            if(count > 1e5){
                std::cerr << "Warning: Max timesteps for dispersion reached"  << std::endl;
                break;
            }
        }
        glm::dvec3 com;
        for(i=0;i<NCELLS;i++){
            com = Cells[i].GetCOM();
            for(j=0;j<Cells[i].NV;j++){
                Cells[i].Positions[j] -= com;
                Cells[i].Positions[j] += centers[i];
            }
        }
    }

    void Tissue::RetractingForceUpdate(){
        std::vector<std::thread> threads;
        int i;
        for(i=0;i<NCELLS;i++){
            threads.push_back(std::thread(&DPM3D::Tissue::CellRetractingUpdate,this,i));
        }
        for(i=0;i<NCELLS;i++){
            threads[i].join();
        }
    }

    void Tissue::CellRetractingUpdate(int ci){
        glm::dvec3 comi = Cells[ci].GetCOM(),comj, rij;
        int vi;
        for(int j=0;j<NCELLS;j++){
            if(ci != j){
                comj = Cells[j].GetCOM();
                FindOverlaps(ci,j);
                for(vi=0; vi<Cells[ci].NV;vi++){
                    if(overlaps[vi]){
                        Cells[ci].Forces[vi] += Kc * glm::normalize(rij);
                    }
                }
            }
        }
    }


    void Tissue::UpdateShapeForces(){
        for(int i=0;i<NCELLS;i++){
            Cells[i].ShapeForceUpdate();
        }
    }


    void Tissue::EulerUpdate(int steps, double dt){
        int step;
        for(step=0;step<steps;step++)
        {   
            UpdateShapeForces();
            RetractingForceUpdate();
            EulerUpdate(dt);
        }
    }

    void Tissue::EulerUpdate(double dt){
        for(int i=0;i<NCELLS;i++){
            Cells[i].EulerUpdate(dt);
        }
    }

    void Tissue::FindOverlaps(int ci, int cj){
        int vi;
        glm::dvec3 point, com = Cells[cj].GetCOM();
        overlaps.resize(Cells[ci].NV);
        for(vi = 0; vi < Cells[ci].NV; vi++){
            point = Cells[ci].Positions[vi];
            point -= L * round(point/L);
            point += L * round(com/L);
            overlaps[vi] = Cells[cj].pointInside(point);
        }
    }

}
