#include <cmath>
#include <Tissue.hpp>
#include <vector>
#include <iostream>
#include <glm/glm.hpp>
#include <glm/vec3.hpp>
#include <glm/gtx/norm.hpp>

namespace DPM3D{
    Tissue::Tissue(std::vector<DPM3D::Cell> _cells, double _phi0){
        Cells = _cells;
        phi0 = _phi0;
        NCELLS = Cells.size();
    }
    
    void Tissue::InteractingForceUpdate(){
        int i,j,vi,vj;
        glm::dvec3 comi, comj;
        double avgdiameter,dist,ftmp;
        for(i=0;i<NCELLS;i++){
            comi = Cells[i].GetCOM();
            for(j=0;j<NCELLS;j++){
                if(i!=j){
                    comj = Cells[j].GetCOM();
                    for(vi=0;vi<Cells[i].NV;vi++){
                        for(vj=0;vj<Cells[j].NV;vj++){
                            avgdiameter = distance(comi,Cells[i].Positions[vi]) + distance(comj,Cells[j].Positions[vi]);
                            dist = distance(comi,comj);
                            if(dist<avgdiameter){
                                ftmp = (Kc/2)*pow((dist/avgdiameter),2);
                                Cells[i].Forces[vi] += ftmp*glm::normalize(Cells[j].Positions[vj] - Cells[i].Positions[vi]);
                                Cells[j].Forces[vj] -= ftmp*glm::normalize(Cells[j].Positions[vj] - Cells[i].Positions[vi]);
                            }
                        }
                    }
                }
            }
        }
    }

}