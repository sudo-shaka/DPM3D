/*
 * =====================================================================================
 *
 *       Filename:  ECM.cpp
 *
 *    Description:
 *
 *        Version:  1.0
 *        Created:  09/01/2022 02:13:14 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (),
 *   Organization:
 *
 * =====================================================================================
 */
#include "ECM.hpp"
#include <Cell3D.hpp>
#include <glm/glm.hpp>
#include <glm/vec3.hpp>
#include <iostream>

namespace DPM3D{
  ECM::ECM(std::vector<glm::dvec3> Points){
    NP = Points.size();
    Forces.resize(NP);
    Positions.resize(NP);

    for(int i = 0; i < NP; i++){
      Forces[i] = {0,0,0};
      Positions[i] = Points[i];
    }
  }

  ECM::ECM(double Z, double L, int _NP){
    NP = sqrt(_NP);
    double len = sqrt(L);
    double dl = len/NP;
    int i,j;
    double x,y=0.0;
    glm::dvec3 tmp;

    for(i=0;i<NP;i++){
      x = 0.0;
      y += dl;
      for(j=0;j<NP;j++){
        x += dl;
        tmp = {x,y,Z};
        Positions.push_back(tmp);
        tmp = {0,0,0};
        Forces.push_back(tmp);
      }
    }

    NP = Positions.size();
    std::cout << NP << std::endl;
  }

  ECM::ECM(){
    Positions.resize(0);
    Forces.resize(0);
    NP = 0;
  }

  void ECM::StretchX(double scale){
    for(int i=0; i < NP; i++){
      Positions[i].x *= 1.0+scale;
    }
  }
  void ECM::StretchY(double scale){
    for(int i=0; i < NP; i++){
      Positions[i].y *= 1.0+scale;
    }
  }
  void ECM::StretchZ(double scale){
    for(int i=0; i < NP; i++){
      Positions[i].z *= 1.0+scale;
    }
  }
}
