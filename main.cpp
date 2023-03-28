#include <iostream>
#include "OpenMesh/Core/Mesh/Handles.hh"
#include "src/common.h"
#include "Eigen/Core"
#include "Eigen/Dense"
#include "src/FrameField.h"
#include "OpenMesh/Core/IO/MeshIO.hh"
#include "string"

#define LEN 1.0

using namespace std;

void generateFrameFieldViz(Mesh& mesh, FrameField& f, string path){
    Mesh m;
    for(auto it = mesh.faces_sbegin();it!=mesh.faces_end();it++){
        auto fh = *it;
        auto cf = f.getCrossFields(fh);
        Eigen::Vector3d c = mesh.calc_centroid(fh);
        auto cvh = m.add_vertex_dirty(c);
        for(auto d:cf){
            auto vh1 = m.add_vertex_dirty(c+LEN*d);
            auto vh2 = m.add_vertex_dirty(c+LEN*d);
            m.add_face(vector<OpenMesh::VertexHandle>{vh1,vh2,cvh});
        }
    }
    OpenMesh::IO::write_mesh(m,path);
}

int main(int, char**) {

    Mesh m;
    OpenMesh::IO::read_mesh(m,"../data/woody.obj");
    FrameField f(m);
    f.generateBoundaryConstrain();
    f.generateFrameField();
    // 生成可视化
    generateFrameFieldViz(m, f, "../data/cf.obj");

    // 奇异点判断
    for(auto vh:m.vertices()){
        auto s = f.isSingular(vh);
        if(s > 0 && s != 4){
            cout << vh << endl;
        }
    }
}
