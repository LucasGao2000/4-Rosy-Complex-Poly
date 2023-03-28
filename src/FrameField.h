#pragma once

#include "OpenMesh/Core/Mesh/Handles.hh"
#include "vector"
#include "common.h"
#include "complex"
#include "unordered_map"

/**
* 标架场类
* 根据OpenMesh的
*/
class FrameField{
public:
    FrameField(Mesh& _mesh);

    // 根据边界加入约束
    void generateBoundaryConstrain();

    // 生成标架场
    void generateFrameField();

    // 获取标架场
    std::vector<Eigen::Vector3d> getCrossFields(OpenMesh::FaceHandle fh);
    std::vector<Eigen::Vector3d> getCrossFields(int _fid);

    // 是否是奇异点
    int isSingular(OpenMesh::VertexHandle vh);

private:
    // Openmesh 网格数据结构
    Mesh* mesh;

    // 局部坐标系
    std::vector<std::vector<Eigen::Vector3d>> localCoords;

    // 所有的面上的complex
    std::vector<std::complex<double>> allComplex;

    // 约束
    std::unordered_map<OpenMesh::FaceHandle, std::complex<double>> constrains;

    // 生成局部坐标系
    void generateLocalCoordinates();

};