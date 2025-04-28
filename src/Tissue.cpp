#define GLM_ENABLE_EXPERIMENTAL
#include <cassert>
#include <cmath>
#include "Tissue.hpp"
#include <vector>
#include <array>
#include <string>
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
            volume+=(4.0f/3.0f)*M_PI*pow(Cells[i].r0,3);
        }
        L = cbrt(volume)/phi0;
        PBC = true;
        attactionMethod.assign("General");
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
                ri = Cells[i].r0*2;
                for(j=0;j<NCELLS;j++){
                    if(j != i){
                        xj = X[j];
                        yj = Y[j];
                        rj = Cells[j].r0*2;
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

    void Tissue::InteractingUpdate(){
        std::vector<std::thread> threads;
        int i;
        for(i=0;i<NCELLS;i++){
            threads.push_back(std::thread(&DPM3D::Tissue::CellInteractingUpdate,this,i));
        }
        for(auto& th : threads){
            th.join();
        }
    }


    void Tissue::UpdateJunctions(){
      std::vector<std::thread> threads;
      for(int i=0; i<NCELLS;i++){
        threads.push_back(std::thread(&DPM3D::Cell::FindJunctions,&this->Cells[i]));
      }
      for(auto& th : threads){
        th.join();
      }
    }

    void Tissue::JunctionSlipForceUpdate(int ci){
      UpdateJunctions();
      double dist, l0 = sqrt((4*Cells[ci].a0)/sqrt(3));
      for(int vi = 0;vi < Cells[ci].NV; vi++){
        if(Cells[ci].isJunction[vi]){
          for(int cj = 0; cj < NCELLS; cj++){
            if(ci != cj){
              for(int vj=0;vj<Cells[cj].NV;vj++){
                if(Cells[cj].isJunction[vj]){
                  glm::dvec3 rij = Cells[cj].Positions[vj] - Cells[ci].Positions[vi];
                  if(PBC){rij-= L*round(rij/L);}
                  dist = sqrt(dot(rij,rij));
                  if(dist < l0){
                      Cells[ci].Forces[vi] -= Kat * 0.5 * ((dist/l0) - 1.0) * (rij/dist);
                  }
                }
              }
            }
          }
        }
      }
    }

    void Tissue::JunctionCatchForceUpdate(int ci){
      UpdateJunctions();
      for(int vi=0;vi<Cells[ci].NV;ci++){
        glm::dvec3 *minrij = NULL;
        double mindist=Cells[ci].r0*2;
        if(!Cells[ci].isJunction[vi]){
          continue;
        }
        for(int cj=0;cj<NCELLS;cj++){
          if(ci == cj){
            continue;
          }
          for(int vj=0;vj<Cells[cj].NV;vj++){
            if(!Cells[cj].isJunction[vj]){
              continue;
            }
            glm::dvec3 rij = Cells[cj].Positions[vj] - Cells[ci].Positions[vi];
            if(PBC){rij -= L*round(rij/L);};
            double dist = sqrt(dot(rij,rij));
            if(dist < mindist){
              mindist = dist;
              minrij = &rij;
            }
          }
          if(minrij){
            glm::dvec3 delta = *minrij;
            //Cells[ci].Forces[vi] += delta * Kat * 0.5 * (pow(mindist,2));
            if(mindist != 0)
              Cells[ci].Forces[vi] += glm::normalize(delta) * (Kat * 0.5)/(pow(mindist,2));
          }
        }
      }
    }

    void Tissue::GeneralAttraction(int ci){
      double dist;
      double l0 = sqrt((4*Cells[ci].a0)/sqrt(3));
//      glm::dvec3 com = Cells[ci].GetCOM();
      for(int vi=0; vi < Cells[ci].NV; vi++){
        for(int cj=0;cj<NCELLS;cj++){
          if(ci!=cj){
            for(int vj=0;vj<Cells[cj].NV;vj++){
              glm::dvec3 rij = Cells[cj].Positions[vj] - Cells[ci].Positions[vi];
              if(PBC){
                rij -= L*round(rij/L);
              }
              dist = sqrt(glm::dot(rij,rij));
              //if(dist < l0 && distance(com,Cells[cj].Positions[vj]) > distance(com,Cells[ci].Positions[vi])){
              if(dist < l0){
                Cells[ci].Forces[vi] -= Kat * 0.5 * ((dist/l0) - 1.0) *
                  glm::normalize(Cells[cj].Positions[vj]-Cells[ci].Positions[vi]);
              }
            }
          }
        }
       }
    }

    void Tissue::CellInteractingUpdate(int ci){
      std::vector<int> tri{0,0,0}, trij{0,0,0};
      int fi,cj,fj,NT = Cells[ci].ntriangles;
      glm::dmat3 PI, PJ;
      glm::dvec3 rij, normali, normalj, FaceCenterI, FaceCenterJ;
      double l0 = sqrt((4*Cells[ci].a0)/sqrt(3));
      double dist;
      glm::dvec3 com = Cells[ci].GetCOM();

      if(Kat != 0){
        if(attactionMethod == "JunctionSlip"){
          JunctionSlipForceUpdate(ci);
        }
        else if(attactionMethod == "JunctionCatch"){
          JunctionCatchForceUpdate(ci);
        }
        else if(attactionMethod == "General"){
          GeneralAttraction(ci);
        }
        else{
          std::cerr << attactionMethod << " is not a valid attaction method\n" <<
            "please use JunctionSlip, JunctionCatch, or General" << std::endl;
          exit(1);
        }
      }

      //Repulsive Faces
      for(fi=0;fi<NT;fi++){
        tri = Cells[ci].FaceIndices[fi];
        PI[0] = Cells[ci].Positions[tri[0]];
        PI[1] = Cells[ci].Positions[tri[1]];
        PI[2] = Cells[ci].Positions[tri[2]];
        normali = glm::cross((PI[1]-PI[0]),(PI[2]-PI[0]));
        FaceCenterI = (PI[0]+PI[1]+PI[2])/3.0;
        for(cj=0;cj<NCELLS;cj++){
          if(cj != ci){
            for(fj=0;fj<Cells[cj].ntriangles;fj++){
              trij = Cells[cj].FaceIndices[fj];
              PJ[0] = Cells[cj].Positions[trij[0]];
              PJ[1] = Cells[cj].Positions[trij[1]];
              PJ[2] = Cells[cj].Positions[trij[2]];
              FaceCenterJ = (PJ[0]+PJ[1]+PJ[2])/3.0;
              rij = FaceCenterJ - FaceCenterI;
              if(PBC){
                rij -= L*round(rij/L);
              }
              dist = sqrt(glm::dot(rij,rij));
              normalj = glm::cross((PJ[1]-PJ[0]),(PJ[2]-PJ[0]));
              if(glm::dot(normali,rij) < 0.0 && dist < l0){
                Cells[ci].Forces[tri[0]] += dist*0.5*Kre*glm::normalize(com-Cells[ci].Positions[tri[0]]);
                Cells[ci].Forces[tri[1]] += dist*0.5*Kre*glm::normalize(com-Cells[ci].Positions[tri[1]]);
                Cells[ci].Forces[tri[2]] += dist*0.5*Kre*glm::normalize(com-Cells[ci].Positions[tri[2]]);
              }
            }
          }
        }
      }
    }

    void Tissue::UpdateShapeForces(){
      std::vector<std::thread> threads;
        for(int i=0;i<NCELLS;i++){
          threads.push_back(std::thread(&DPM3D::Cell::ShapeForceUpdate,&this->Cells[i]));
        }
        for(auto& th : threads){
          th.join();
        }
    }

    void Tissue::StickToSurface(double z, double mindist){
      std::vector<std::thread> threads;
      int i;
      for(i=0;i<NCELLS;i++){
        threads.push_back(std::thread(&DPM3D::Tissue::CellStickToSurface,this,i,z,mindist));
      }
      for(auto& th : threads){
        th.join();
      }
    }

    void Tissue::CellStickToSurface(int ci, double z, double mindist){
      Cells[ci].StickToSurface(z,mindist);
    }

    void Tissue::EulerUpdate(int steps, double dt){
        int step;
        for(step=0;step<steps;step++)
        {
            UpdateShapeForces();
            InteractingUpdate();
            EulerUpdate(dt);
        }
    }

    void Tissue::applyShearStress(std::array<double,3> _fluidVelocity, double viscosity, double velocityGradient, double stiffness){
      glm::dvec3 fluidVelocity = {_fluidVelocity[0],_fluidVelocity[1],_fluidVelocity[2]};
      if(stiffness == 0)
        return;
      if(fluidVelocity.x == 0 && fluidVelocity.y == 0 && fluidVelocity.z ==0)
        return;
      if(viscosity==0)
        return;

      for(auto cell : Cells){
        int flen = cell.FaceIndices.size();
        for(int fi=0;fi < flen; fi++){
          glm::ivec3 face = {cell.FaceIndices[fi][0],cell.FaceIndices[fi][1],cell.FaceIndices[fi][2]};
          glm::dvec3 p1 = cell.Positions[face[0]];
          glm::dvec3 p2 = cell.Positions[face[1]];
          glm::dvec3 p3 = cell.Positions[face[2]];
          glm::dvec3 normal = glm::normalize(glm::cross((p2-p1),(p3-p1)));
          glm::dvec3 velocityParrallel = fluidVelocity - glm::dot(fluidVelocity,normal) * normal;
          glm::dvec3 shearStress = viscosity * velocityGradient * velocityParrallel;
          double area = cell.GetArea(fi);
          glm::dvec3 shearForce = shearStress * area;
          cell.Forces[face[0]] += shearForce * stiffness;
          cell.Forces[face[1]] += shearForce * stiffness;
          cell.Forces[face[2]] += shearForce * stiffness;
        }
      }
    }

    void Tissue::InteractECM(DPM3D::ECM ECM){
      int pi,ci,vi;
      double l0,dist;
      glm::dvec3 rij,COM;
      for(pi=0;pi<ECM.NP;pi++){
        ECM.Forces[pi] *= 0.0;
        for(ci=0;ci<NCELLS;ci++){
          COM = Cells[ci].GetCOM();
          l0 = sqrt((4*Cells[ci].a0)/sqrt(3));
          for(vi=0;vi<Cells[ci].NV;vi++){
            dist = distance(Cells[ci].Positions[vi],ECM.Positions[pi]);
            if(dist < l0 && distance(ECM.Positions[pi],COM) > distance(Cells[ci].Positions[vi],COM)){
              rij = COM - ECM.Positions[pi];
              ECM.Forces[pi] -= Cells[ci].Ks * ((1-dist/(l0))/l0) *glm::normalize(rij);
              Cells[ci].Positions[vi] += Cells[ci].Ks * ((1-dist/(l0))/l0) *glm::normalize(rij);
            }
          }
        }
      }
    }

    void Tissue::EulerUpdate(double dt){
        for(int i=0;i<NCELLS;i++){
            Cells[i].EulerUpdate(dt);
        }
    }

    std::vector<bool> Tissue::FindOverlaps(int ci, int cj){
        std::vector<bool> overlaps;
        glm::dvec3 point, com = Cells[cj].GetCOM();
        overlaps.resize(Cells[ci].NV);
        for(int vi = 0; vi < Cells[ci].NV; vi++){
            point = Cells[ci].Positions[vi];
            if(PBC){
              point -= L * round(point/L);
              point += L * round(com/L);
            }
            overlaps[vi] = Cells[cj].pointInside(point);
        }
        return overlaps;
    }

    std::vector<std::vector<double>> Tissue::GetVesselPosition(int ci){
      std::vector<std::vector<double>> vesselPosition;
      if(!PBC){
        return vesselPosition;
      }
      DPM3D::Cell c = Cells[ci];
      vesselPosition.resize(3);
      for(int i=0; i<3;i++){
        vesselPosition[i].resize(c.NV);
      }
      double scale = (2.0*M_PI)/L, radius = L/(2*M_PI);
      for(int vi=0;vi<c.NV;vi++){
        double theta = c.Positions[vi].x * scale;
        vesselPosition[0][vi] = (radius-c.Positions[vi].z) * cos(theta);
        vesselPosition[2][vi] = (radius-c.Positions[vi].z) * sin(theta);
        vesselPosition[1][vi] = c.Positions[vi].y;
      }
      return vesselPosition;
    }
}
