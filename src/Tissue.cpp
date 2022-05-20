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
            volume+=Cells[i].v0;
        }
        L = pow(volume,(1/3))/phi0;
    }

    void Tissue::RetractingForceUpdate(){
        int i,j,vi;//vj;
        glm::dvec3 comi, comj;
        double avgdiameter,dx,dy,dz,dist,ftmp;
        for(i=0;i<NCELLS;i++){
            comi = Cells[i].GetCOM();
            for(j=0;j<NCELLS;j++){
                comj = Cells[j].GetCOM();
                if(i!=j){
                    FindOverlaps(i,j);
                    for(vi=0;vi<Cells[i].NV;vi++){
                        /*dx = Cells[i].Positions[vi].x-comi.x;
                        dy = Cells[i].Positions[vi].y-comi.y;
                        dz = Cells[i].Positions[vi].z-comi.z;
                        dx -= L*round(dx/L);
                        dy -= L*round(dy/L);
                        dz -= L*round(dz/L);
                        dist = sqrt(dx*dx + dy*dy + dz*dz);*/
                        if(overlaps[vi]){
                            dx = Cells[i].Positions[vi].x-comi.x;
                            dy = Cells[i].Positions[vi].y-comi.y;
                            dz = Cells[i].Positions[vi].z-comi.z;
                            dx -= L*round(dx/L);
                            dy -= L*round(dy/L);
                            dz -= L*round(dz/L);
                            dist = sqrt(dx*dx + dy*dy + dz*dz);
                            ftmp = Kc*(1-dist)/Cells[i].r0;
                            Cells[i].Forces[vi] -= ftmp * glm::normalize(comi - Cells[i].Positions[vi]);
                        }
                        /*else if(dist < 1.0){
                            ftmp = Kc*(1-dist)/Cells[i].r0;
                            Cells[i].Forces[vi] += ftmp * glm::normalize(comj - Cells[i].Positions[vi]);
                        }*/
                    }
                }
            }
        }
    }


    void Tissue::EulerUpdate(int steps, double dt){
        int i,step;
        for(step=0;step<steps;step++)
        {
            for(i=0;i<NCELLS;i++){
                Cells[i].ShapeForceUpdate();
            }
            RetractingForceUpdate();
            for(i=0;i<NCELLS;i++){
                Cells[i].EulerUpdate(dt);
            }
        }
    }

    void Tissue::FindOverlaps(int ci, int cj){
        int vi;
        //double dx,dy,dz;
        //glm::dvec3 point, com = Cells[cj].GetCOM();
        overlaps.resize(Cells[ci].NV);
        for(vi = 0; vi < Cells[ci].NV; vi++){
            /*dx = Cells[i].Positions[vi].x - com.x;
            dy = Cells[i].Positions[vi].y - com.y;
            dz = Cells[i].Positions[vi].y - com.y;
            dx = L*round(dx/L);
            dy = L*round(dy/L);
            dz = L*round(dz/L);
            point.x = dx;
            point.y = dy;
            point.z = dz;*/
            overlaps[vi] = Cells[cj].pointInside(Cells[ci].Positions[vi]);
        }
    }

}
