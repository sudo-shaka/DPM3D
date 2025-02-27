#ifndef __ECM__
#define __ECM__

#include<glm/vec3.hpp>
#include<vector>

namespace DPM3D{
  struct ECM{
      int NP;
      std::vector<glm::dvec3> Positions;
      std::vector<glm::dvec3> Forces;

      ECM(std::vector<glm::dvec3> Points);
      ECM(double Z, double L, int NP);
      ECM();

      void StretchX(double scale);
      void StretchY(double scale);
      void StretchZ(double scale);

  };
}
#endif
