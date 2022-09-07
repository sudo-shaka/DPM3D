/*
 * =====================================================================================
 *
 *       Filename:  ECM.hpp
 *
 *    Description: Extracellular Matrix
 *
 *        Version:  1.0
 *        Created:  09/01/2022 02:18:32 PM
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  YOUR NAME (),
 *   Organization:
 *
 * =====================================================================================
 */

#include<vector>
#include<glm/vec3.hpp>

namespace DPM3D{
  struct ECM{
      int NP;
      std::vector<glm::dvec3> Positions;
      std::vector<glm::dvec3> Forces;

      ECM(std::vector<glm::dvec3> Points);
      ECM(double Z, double L, int NP);
  };
}
