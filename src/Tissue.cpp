#define GLM_ENABLE_EXPERIMENTAL
#include <algorithm>
#include <glm/common.hpp>
#include <glm/vector_relational.hpp>
#include "Tissue.hpp"
#include <glm/geometric.hpp>
#include <vector>
#include <array>
#include <string>
#include <iostream>
#include <glm/glm.hpp>
#include <glm/vec3.hpp>
#include <glm/gtx/norm.hpp>
#include <thread>
#include <fstream>

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
        double ri,rj,yi,yj,xi,xj,dx,dy,dist, U;
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
                        dx -= L*std::round(dx/L);
                        dy = yj-yi;
                        dy -= L*std::round(dy/L);
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
                        rij -= L*glm::round(rij/L);
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
      double dist, l0 = Cells[ci].l0;
      for(int vi = 0;vi < Cells[ci].NV; vi++){
        if(Cells[ci].isJunction[vi]){
          for(int cj = 0; cj < NCELLS; cj++){
            if(ci != cj){
              for(int vj=0;vj<Cells[cj].NV;vj++){
                if(Cells[cj].isJunction[vj]){
                  glm::dvec3 rij = Cells[cj].Positions[vj] - Cells[ci].Positions[vi];
                  if(PBC){rij-= L*glm::round(rij/L);}
                  dist = sqrt(dot(rij,rij));
                  if(dist < l0){
                    double ftmp = dist/l0 * Kat;
                    double lifetime = std::exp(std::fabs(ftmp)/f0);
                    ftmp /= lifetime;
                    Cells[ci].Forces[vi] += 0.5 * ftmp * glm::normalize(Cells[cj].Positions[vj] - Cells[ci].Positions[vi]);
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
      double dist, l0 = Cells[ci].l0;
      for(int vi = 0;vi < Cells[ci].NV; vi++){
        if(Cells[ci].isJunction[vi]){
          for(int cj = 0; cj < NCELLS; cj++){
            if(ci != cj){
              for(int vj=0;vj<Cells[cj].NV;vj++){
                if(Cells[cj].isJunction[vj]){
                  glm::dvec3 rij = Cells[cj].Positions[vj] - Cells[ci].Positions[vi];
                  if(PBC){rij-= L*glm::round(rij/L);}
                  dist = sqrt(dot(rij,rij));
                  if(dist < l0){
                    double ftmp = dist/l0 * Kat;
                    double lifetime = std::exp(-std::fabs(ftmp)/f0) * 0.5 * std::exp(-std::pow((std::fabs(ftmp)-f0)/f0,2));
                    if(lifetime < 1e-8) continue;
                    ftmp /= lifetime;
                    Cells[ci].Forces[vi] += 0.5 * ftmp * glm::normalize(Cells[cj].Positions[vj] - Cells[ci].Positions[vi]);
                  }
                }
              }
            }
          }
        }
      }
    }

    void Tissue::GeneralAttraction(int ci){
      double dist;
      double l0 = Cells[ci].l0;
      for(int vi=0; vi < Cells[ci].NV; vi++){
        for(int cj=0;cj<NCELLS;cj++){
          if(ci==cj) continue;
          for(int vj=0;vj<Cells[cj].NV;vj++){
            glm::dvec3 rij = Cells[cj].Positions[vj] - Cells[ci].Positions[vi];
            if(PBC){
              rij -= L*glm::round(rij/L);
            }
            dist = sqrt(glm::dot(rij,rij));
            if(dist > l0) continue;
            Cells[ci].Forces[vi] -= Kat * 0.5 * ((dist/l0) - 1.0) *
              glm::normalize(Cells[cj].Positions[vj]-Cells[ci].Positions[vi]);
          }
        }
      }
    }

    void Tissue::CellInteractingUpdate(int ci){
      std::vector<int> tri{0,0,0}, trij{0,0,0};
      if(Kat != 0){
        if(attactionMethod == "JunctionSlip") JunctionSlipForceUpdate(ci);
        else if(attactionMethod == "JunctionCatch") JunctionCatchForceUpdate(ci);
        else if(attactionMethod == "General") GeneralAttraction(ci);
        else if(attactionMethod == "PolarizedJunction") PolarizedJunctionForces(ci);
        else std::cerr << attactionMethod
              << " is not a valid attaction method\n" 
              << "Please use JunctionSlip, JunctionCatch, or General" 
              << "... no attraction force update performed."
              << std::endl;
      }
       //prevent overlaps
      if(Kre == 0) return;

      for(int cj=0; cj < NCELLS;cj++){
        if(cj == ci) continue;
        std::vector<double> windingNumbers = findWindingNumber(ci,cj);
        for(int vi=0;vi<Cells[ci].NV;vi++){
          Cells[ci].Forces[vi] += std::fabs(windingNumbers[vi]) * 0.5 * Kre * glm::normalize(Cells[ci].COM - Cells[ci].Positions[vi]);
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

    void Tissue::StickToSurfaceStretch(double z, double mindist, double xx, double xy){
      std::vector<std::thread> threads;
      for(int i=0;i<NCELLS;i++){
        threads.push_back(std::thread(&DPM3D::Cell::StickToSurfaceStretch,&this->Cells[i],z,mindist,xx,xy));
      }
      for(auto& thread : threads){
        thread.join();
      }
    }

    void Tissue::PolarizedJunctionForces(int ci){
      Cells[ci].FindJunctions();
      std::vector<double> x,y;
      std::vector<std::vector<glm::dvec3>> rijs;
      for(int vi=0;vi<Cells[ci].NV;vi++){
        if(Cells[ci].isJunction[vi]){
          glm::dvec3 vert=Cells[ci].Positions[vi];
          x.push_back(vert.x);
          y.push_back(vert.y);
          std::vector<glm::dvec3> riji;
          for(auto& cell : Cells){
            for(int vj=0;vj<cell.NV;vj++){
              if(distance(cell.Positions[vj],vert) < Cells[ci].l0){
                glm::dvec3 rij = cell.Positions[vj] - Cells[ci].Positions[vi];
                riji.push_back(rij);
              }
            }
          }
          rijs.push_back(riji);
        }
      }
      double minx = *std::min_element(x.begin(),x.end());
      double maxx = *std::max_element(x.begin(),x.end());
      double miny = *std::min_element(y.begin(),y.end());
      double maxy = *std::min_element(y.begin(),y.end());
      auto slipForce = [](const glm::dvec3 rij,const double l0,const double f){
        double ftmp = std::sqrt(glm::dot(rij,rij))/l0;
        double lifetime = std::exp(std::fabs(ftmp)/f);
        if(lifetime < 1e-8) return glm::dvec3(0.0);
        ftmp /= lifetime;
        return 0.5*std::abs(ftmp)*glm::normalize(rij);
      };
      auto catchForce = [](const glm::dvec3 rij,const double l0, const double f){
        double ftmp = std::sqrt(glm::dot(rij,rij)) / l0;
        double lifetime = std::exp(-std::fabs(ftmp)/f) * 0.5 * std::exp(-std::pow((std::fabs(ftmp)-f)/f,2));
        if(lifetime < 1e-8) return glm::dvec3(0.0);
        ftmp /= lifetime;
        return 0.5*std::abs(ftmp)*glm::normalize(rij);
      };
      
      int ji = 0;
      for(int vi=0;vi<Cells[ci].NV;vi++){
        if(Cells[ci].isJunction[vi]){
          glm::dvec3 vert=Cells[ci].Positions[vi];
          std::vector<glm::dvec3> R = rijs[ji];
          if(vert.x < (Cells[ci].COM.x + minx)/2.0 || vert.x > (Cells[ci].COM.x + maxx)/2.0){
            for(int ri=0;ri<(int)R.size();ri++){
              auto force = Kat * catchForce(R[ri],Cells[ci].l0, f0);
              if(!glm::all(glm::isnan(force))) Cells[ci].Forces[vi] += force;
            }
          }
          if(vert.y < (Cells[ci].COM.y + miny)/2.0 || vert.y > (Cells[ci].COM.y+maxy)/2.0){
            for(int ri=0;ri<(int)R.size();ri++){
              auto force = Kat * slipForce(R[ri],Cells[ci].l0, f0);
              if(!glm::all(glm::isnan(force))) Cells[ci].Forces[vi] += force;
            }
          }
          ji++;
        }
      }
    }

    void Tissue::PolarizedJunctionUpdate(){
      std::vector<std::thread> threads;
      for(int i=0;i<NCELLS;i++){
        threads.push_back(std::thread(&DPM3D::Tissue::PolarizedJunctionForces,this,i));
      }
      for(auto& th : threads){
        th.join();
      }
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
          COM = Cells[ci].COM;
          l0 = Cells[ci].l0;
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

    std::vector<double> Tissue::findWindingNumber(int ci, int cj){
      std::vector<double> windingNumber(Cells[ci].NV,0.0);
      if(ci == cj){
        std::cerr << "Error: Trying to find if cell is inside itself" << std::endl;
        return windingNumber;
      }
      
      glm::dvec3 shift(0.0); 
      if(PBC){
        shift = L * glm::round((Cells[ci].COM-Cells[cj].COM)/L);
      }
      //AABB check first
      glm::dvec3 minA(FLT_MAX), maxA(-FLT_MAX), minB(FLT_MAX), maxB(-FLT_MAX);
      for(const auto& v : Cells[ci].Positions){
        minA = glm::min(minA,v);
        maxA = glm::max(maxA,v);
      }

      for(const auto& v : Cells[cj].Positions){
        minB = glm::min(minB,v+shift);
        maxB = glm::max(maxB,v+shift);
      }

      bool overlapX = maxA.x >= minB.x && minA.x <= maxB.x;
      bool overlapY = maxA.y >= minB.y && minA.y <= maxB.y;
      bool overlapZ = maxA.z >= minB.z && minA.z <= maxB.z;
      bool BBoverlap = overlapX && overlapY && overlapZ;
      if(!BBoverlap) return windingNumber;
      for(int vi=0;vi<Cells[ci].NV;vi++){
        glm::dvec3 point = Cells[ci].Positions[vi] - shift;
        windingNumber[vi] = Cells[cj].WindingNumberOf(point);
      }
      return windingNumber;
    }

    std::vector<bool> Tissue::FindOverlaps(int ci, int cj){
      std::vector<bool> overlaps(Cells[ci].NV);
      std::vector<double> windingNumber = findWindingNumber(ci,cj);
      for(int vi = 0; vi < Cells[ci].NV; vi++){
        overlaps[vi] = std::fabs(windingNumber[vi]) > 0.9;
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

    void Tissue::ExportPositions(std::string filename){
      std::ofstream file;
      file.open(filename);
      for(int i=0; i<NCELLS;i++){
        file << "Cell_" << i << ",X,";
        for(auto& vert : Cells[i].Positions){
          file << vert.x << ",";
        }
        file << "\nCell_" << i << ",Y,";
        for(auto& vert : Cells[i].Positions){
          file << vert.y << ",";
        }
        file << "\nCell_" << i << ",Z,";
        for(auto& vert : Cells[i].Positions){
          file << vert.z << ",";
        }
        file << std::endl;
      }
      file.close();
    }

    void Tissue::ExportForces(std::string filename){
      std::ofstream file;
      file.open(filename);
      for(int i=0; i<NCELLS;i++){
        file << "Cell_" << i << ",X,";
        for(auto& force : Cells[i].Forces){
          file << force.x << ",";
        }
        file << "\nCell_" << i << ",Y,";
        for(auto& force : Cells[i].Forces){
          file << force.y << ",";
        }
        file << "\nCell_" << i << ",Z,";
        for(auto& force : Cells[i].Forces){
          file << force.z << ",";
        }
        file << std::endl;
      }
      file.close();
    }
}
