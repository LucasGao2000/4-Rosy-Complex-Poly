#include "FrameField.h"
#include "Eigen/src/Core/Matrix.h"
#include "OpenMesh/Core/Mesh/Handles.hh"
#include <cassert>
#include <complex>
#include <vector>
#include "assert.h"
#include "Eigen/Eigen"
#include "cmath"

typedef Eigen::SparseMatrix<std::complex<double>> SpMat;
typedef Eigen::Triplet<std::complex<double>> T;

using namespace std;

FrameField::FrameField(Mesh& _mesh){
    this->mesh = &_mesh;
    allComplex.resize(mesh->n_faces(), complex<double>(0,0));
    localCoords.resize(mesh->n_faces(), vector<Eigen::Vector3d>(2));
    if(!mesh->has_face_status()) mesh->request_face_status();
    if(!mesh->has_edge_status()) mesh->request_edge_status();
    if(!mesh->has_vertex_status()) mesh->request_vertex_status();
    if(!mesh->has_face_normals()) mesh->request_face_normals();
    // 生成局部坐标
    generateLocalCoordinates();
}

void FrameField::generateLocalCoordinates(){
    if(!mesh) return;

    for(auto it = mesh->faces_sbegin();it!=mesh->faces_end();it++){
        auto fh = *it;
        vector<OpenMesh::VertexHandle> vhs;
        for(auto vh:fh.vertices_ccw()){
            vhs.push_back(vh);
        }
        Eigen::Vector3d p1 = mesh->point(vhs[0]);
		Eigen::Vector3d p2 = mesh->point(vhs[1]);
		Eigen::Vector3d p3 = mesh->point(vhs[2]);

		Eigen::Vector3d X = (p2 - p1).normalized();
		Eigen::Vector3d n = X.cross((p3 - p1).normalized());
		Eigen::Vector3d Y = n.cross(X);
		n.normalize();
		Y.normalize();
        // 保存局部坐标的两个坐标系
        localCoords[fh.idx()] = {X,Y}; 
    }
}

void FrameField::generateBoundaryConstrain(){
    // 遍历所有边界面
    for(auto it = mesh->edges_sbegin(); it!=mesh->edges_end();it++){
        if(!it->is_boundary()) continue;
        auto eh = *it;
        OpenMesh::SmartFaceHandle fh;
        if(eh.h0().face().is_valid()) fh = eh.h0().face();
        else fh = eh.h1().face();
        assert(fh.is_valid());
        Eigen::Vector3d N = mesh->calc_face_normal(fh);
        Eigen::Vector3d p0 = mesh->point(eh.v0());
        Eigen::Vector3d p1 = mesh->point(eh.v1());
        Eigen::Vector3d dir = (p0-p1).normalized();
        // 垂直于边界的向量
        Eigen::Vector3d pp = dir.cross(N);
        // 设置约束
        double cos = pp.dot(localCoords[fh.idx()][0]);
        double sin = pp.dot(localCoords[fh.idx()][1]);
       
        // cos + i sin
        constrains[fh] = std::complex<double>(cos, sin);
        allComplex[fh.idx()] = std::complex<double>(cos, sin);
    }
}


void FrameField::generateFrameField(){
    if(!mesh || mesh->n_faces() == 0) return;

    vector<T> triples;
    int cnt = 0;
    for(auto it = mesh->faces_sbegin(); it!=mesh->faces_end();it++){
        auto fh = *it;
        for(auto heh:fh.halfedges_ccw()){
            if(!heh.opp().face().is_valid()) continue;
            auto nfh = heh.opp().face();
            Eigen::Vector3d e = (mesh->point(heh.from()) - mesh->point(heh.to())).normalized();
            complex<double> e_f(e.dot(localCoords[fh.idx()][0]),e.dot(localCoords[fh.idx()][1]));
            complex<double> e_g(e.dot(localCoords[nfh.idx()][0]),e.dot(localCoords[nfh.idx()][1]));
            e_f = conj(e_f);
            e_g = conj(e_g);
            triples.emplace_back(cnt, fh.idx(), pow(e_f,4));
            triples.emplace_back(cnt, nfh.idx(), -pow(e_g,4));
            cnt++;
        }
    }

    Eigen::VectorXcd b = Eigen::VectorXcd::Zero(cnt + constrains.size());

    // add contrains
    for(auto it:constrains){
        auto fh = it.first;
        triples.emplace_back(cnt, fh.idx(), 1);
        b[cnt] = pow(it.second,4);
        // cout << it.second << endl;
        cnt++;
    }

    SpMat A(cnt, mesh->n_faces());
    A.setFromTriplets(triples.begin(), triples.end());
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<complex<double>>> solver;
    SpMat AT = A.adjoint();
    solver.compute(AT*A);
	if (solver.info() != Eigen::Success) {
		// decomposition failed
	    cerr << "[ComplexPoly]: solver init failed" << endl;
        cerr << solver.info() << endl;
		return;
	}
	Eigen::VectorXcd x = solver.solve(AT*b);

    // 填入结果
    for(auto it = mesh->faces_sbegin(); it!=mesh->faces_end();it++){
        auto fh = *it;
        if(constrains.contains(fh)){
            continue;
        }
        double r = arg(x[fh.idx()])/4;
        allComplex[fh.idx()] = complex<double>(cos(r),sin(r));
    }
}

vector<Eigen::Vector3d> FrameField::getCrossFields(OpenMesh::FaceHandle fh){
    return getCrossFields(fh.idx());
}

vector<Eigen::Vector3d> FrameField::getCrossFields(int fid){
    vector<Eigen::Vector3d> ret;
    if(allComplex.size() <= fid || fid<0) return ret;
    auto c = allComplex[fid];
    auto& localCoord = localCoords[fid];
    // 其中一条向量
    Eigen::Vector3d dir = c.real() * localCoord[0] + c.imag() * localCoord[1];
    auto N = mesh->calc_face_normal(mesh->face_handle(fid));
    // 第一条向量
    ret.push_back(dir);
    for(int i=0;i<3;i++){
        Eigen::Vector3d d = N.cross(dir);
        dir = d;
        ret.push_back(dir);
    }
    return ret;
}

int FrameField::isSingular(OpenMesh::VertexHandle vh){
    if(!vh.is_valid()) return -1;
    if(mesh->is_boundary(vh)) return -1;
    // 逐个获取标架场
    vector<vector<Eigen::Vector3d>> cfs;
    for(auto fh:mesh->vf_ccw_range(vh)){
        auto cf = getCrossFields(fh);
        cfs.push_back(cf);
    }
    Eigen::Vector3d d = cfs[0][0];
    double angle = 0;
    // 逐个匹配
    for(int i=1;i<cfs.size();i++){
        vector<Eigen::Vector3d>& cf = cfs[i];
        auto cf0 = cf[0];
        auto cf1 = cf[1];
        if(abs(cf0.dot(d)) > abs(cf1.dot(d))){
            angle += acos(clamp(abs(cf0.dot(d)),0.,1.));
            d = cf0;
        }else{
            angle += acos(clamp(abs(cf1.dot(d)),0.,1.));
            d = cf1;
        }
    }
    auto cf0 = cfs[0][0];
    auto cf1 = cfs[0][1];

    // 计算夹角
    auto ohehs = mesh->voh_ccw_range(vh).to_vector();
    double angleAround = 0;
    for(int i=0;i<ohehs.size();i++){
        auto heh1 = ohehs[i];
        auto heh2 = ohehs[(i+1) % ohehs.size()];
        auto d1 = (mesh->point(heh1.to()) - mesh->point(heh1.from()));
        auto d2 = (mesh->point(heh2.to()) - mesh->point(heh2.from()));
        angleAround += acos(clamp(fabs(d1.dot(d2)/(d1.norm() * d2.norm())),0.,1.));
    }
     angle += 2 * M_PI - angleAround;

    // 判断奇异点类型
    int index = 0;
    if (fabs(angle - M_PI * 0.5) < 0.2 * M_PI)
        index = -1;
    else if (fabs(angle + M_PI * 0.5) < 0.2 * M_PI)
        index = 1;
    
    return index + 4;
}