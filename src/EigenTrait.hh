//
// Created by Lucas on 2022/10/29.
//

#ifndef QUADGENERATOR_EIGENTRAIT_HH
#define QUADGENERATOR_EIGENTRAIT_HH

#include "OpenMesh/Core/Geometry/EigenVectorT.hh"
#include "OpenMesh/Core/Mesh/Traits.hh"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Geometry>


namespace OpenMesh {
    struct EigenTraits : OpenMesh::DefaultTraits {
        using Point = Eigen::Vector3d;
        using Normal = Eigen::Vector3d;
        using TexCoord2D = Eigen::Vector2d;
    };
}


#endif //QUADGENERATOR_EIGENTRAIT_HH
