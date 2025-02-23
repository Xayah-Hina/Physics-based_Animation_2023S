#include "Task05.h"

#include <SIM/SIM_Engine.h>
#include <SIM/SIM_Object.h>
#include <SIM/SIM_PositionSimple.h>
#include <SIM/SIM_GeometryCopy.h>
#include <SIM/SIM_DopDescription.h>
#include <PRM/PRM_Template.h>
#include <PRM/PRM_Default.h>
#include <PRM/PRM_Name.h>

#include <GU/GU_PrimPoly.h>
#include <Eigen/Dense>

auto vertex_to_elem(
    const Eigen::MatrixXi& elem2vtx,
    size_t num_vtx)
{
    std::vector<unsigned int> vtx2idx, idx2elem;
    vtx2idx.assign(num_vtx + 1, 0);
    for (int i_elem = 0; i_elem < elem2vtx.rows(); i_elem++)
    {
        for (int i_node = 0; i_node < elem2vtx.cols(); i_node++)
        {
            const int i_vtx = elem2vtx(i_elem, i_node);
            vtx2idx[i_vtx + 1] += 1;
        }
    }
    for (unsigned int i_vtx = 0; i_vtx < num_vtx; ++i_vtx)
    {
        vtx2idx[i_vtx + 1] += vtx2idx[i_vtx];
    }
    unsigned int num_idx = vtx2idx[num_vtx];
    idx2elem.resize(num_idx);
    for (int i_elem = 0; i_elem < elem2vtx.rows(); i_elem++)
    {
        for (int i_node = 0; i_node < elem2vtx.cols(); i_node++)
        {
            const int i_vtx = elem2vtx(i_elem, i_node);
            const unsigned int ind1 = vtx2idx[i_vtx];
            idx2elem[ind1] = i_elem;
            vtx2idx[i_vtx] += 1;
        }
    }
    for (int ivtx = static_cast<int>(num_vtx); ivtx >= 1; --ivtx)
    {
        vtx2idx[ivtx] = vtx2idx[ivtx - 1];
    }
    vtx2idx[0] = 0;
    return std::make_pair(vtx2idx, idx2elem);
}

auto generate_mesh_annulus3(
    float r_small,
    float r_large,
    int ndiv_radius,
    int ndiv_theta)
{
    Eigen::Matrix<int, Eigen::Dynamic, 3, Eigen::RowMajor> tri2vtx(ndiv_radius * ndiv_theta * 2, 3);
    Eigen::Matrix<float, Eigen::Dynamic, 3, Eigen::RowMajor> vtx2xyz((ndiv_radius + 1) * ndiv_theta, 3);
    constexpr float pi = 3.1415926535;
    {
        // make coordinates
        float dr = (r_large - r_small) / static_cast<float>(ndiv_radius);
        float dth = 2.f * pi / static_cast<float>(ndiv_theta);
        for (int ir = 0; ir <= ndiv_radius; ir++)
        {
            for (int ith = 0; ith < ndiv_theta; ith++)
            {
                float rad = dr * static_cast<float>(ir) + r_small;
                float theta = static_cast<float>(ith + (ir % 2) * 0.5) * dth;
                vtx2xyz.row(ir * ndiv_theta + ith) = Eigen::Vector3f(rad * cos(theta), 0.f, rad * sin(theta));
            }
        }
    }
    for (int ir = 0; ir < ndiv_radius; ir++)
    {
        for (int ith = 0; ith < ndiv_theta; ith++)
        {
            int i1 = (ir + 0) * ndiv_theta + (ith + 0) % ndiv_theta;
            int i2 = (ir + 0) * ndiv_theta + (ith + 1) % ndiv_theta;
            int i3 = (ir + 1) * ndiv_theta + (ith + 1) % ndiv_theta;
            int i4 = (ir + 1) * ndiv_theta + (ith + 0) % ndiv_theta;
            if (ir % 2 == 1)
            {
                tri2vtx.row((ir * ndiv_theta + ith) * 2 + 0) = Eigen::Vector3i(i3, i1, i2);
                tri2vtx.row((ir * ndiv_theta + ith) * 2 + 1) = Eigen::Vector3i(i4, i1, i3);
            }
            else
            {
                tri2vtx.row((ir * ndiv_theta + ith) * 2 + 0) = Eigen::Vector3i(i4, i2, i3);
                tri2vtx.row((ir * ndiv_theta + ith) * 2 + 1) = Eigen::Vector3i(i4, i1, i2);
            }
        }
    }
    return std::make_pair(tri2vtx, vtx2xyz);
}

auto lines_of_mesh(
    const Eigen::MatrixXi& elem2vtx,
    int num_vtx)
{
    const auto [vtx2idx, idx2elem] = vertex_to_elem(elem2vtx, num_vtx);
    std::vector<int> _line2vtx;
    for (int i_vtx = 0; i_vtx < num_vtx; ++i_vtx)
    {
        std::set<int> set_connected_points;
        for (unsigned int idx = vtx2idx[i_vtx]; idx < vtx2idx[i_vtx + 1]; ++idx)
        {
            const unsigned int ielem0 = idx2elem[idx];
            for (int inode = 0; inode < elem2vtx.cols(); ++inode)
            {
                const int j_vtx = elem2vtx(ielem0, inode);
                if (j_vtx <= i_vtx) continue;
                set_connected_points.insert(j_vtx);
            }
        }
        for (int j_vtx : set_connected_points)
        {
            _line2vtx.push_back(i_vtx);
            _line2vtx.push_back(j_vtx);
        }
    }
    auto map = Eigen::Map<Eigen::Matrix<int, Eigen::Dynamic, 2, Eigen::RowMajor>>(
        _line2vtx.data(), static_cast<int>(_line2vtx.size() / 2), 2);
    return Eigen::Matrix<int, Eigen::Dynamic, 2, Eigen::RowMajor>(map);
}

void wdw_spring_3d(
    float& w,
    Eigen::Vector3f dw[2],
    const Eigen::Vector3f node2xyz[2],
    const float length_ini,
    const float stiffness)
{
    constexpr int num_node = 2;
    const float length = (node2xyz[0] - node2xyz[1]).norm(); // distance between p0 and p1
    const float C = length - length_ini; // the length differences.
    w = 0.5f * stiffness * C * C; // Hooke's law. energy is square of length difference W=1/2*k*C*C

    // write a few lines of code below to compute the gradient of elastic energy of this spring
    // with respect to the positions of the two end points.

    if (length > 1e-6f)
    {
        // 避免除以0
        Eigen::Vector3f diff = node2xyz[0] - node2xyz[1];
        Eigen::Vector3f grad = stiffness * C * (diff / length);
        dw[0] = grad;
        dw[1] = -grad;
    }
    else
    {
        dw[0].setZero();
        dw[1].setZero();
    }
}

float gradient_descent_energy_minimization(
    Eigen::Matrix<float, Eigen::Dynamic, 3, Eigen::RowMajor>& vtx2xyz,
    const Eigen::Matrix<float, Eigen::Dynamic, 3, Eigen::RowMajor>& vtx2xyz_ini,
    const Eigen::Matrix<int, Eigen::Dynamic, 2, Eigen::RowMajor>& line2vtx,
    float stiffness,
    float mass_point,
    const Eigen::Vector3f& gravity,
    const Eigen::Matrix<float, Eigen::Dynamic, 3, Eigen::RowMajor>& aBCFlag,
    float learning_rate)
{
    // simulation
    const unsigned int num_vtx = vtx2xyz.rows(); // number of vertices
    float W = 0.0; // energy of the system
    Eigen::MatrixX3f gradW = Eigen::MatrixX3f::Zero(num_vtx, 3); // gradient of the energy
    for (int i_line = 0; i_line < line2vtx.rows(); ++i_line)
    {
        // loop over springs
        const int i_vtx0 = line2vtx(i_line, 0); // index of vertex 0
        const int i_vtx1 = line2vtx(i_line, 1); // index of vertex 1
        const float length_ini = (vtx2xyz_ini.row(i_vtx0) - vtx2xyz_ini.row(i_vtx1)).norm(); // initial length
        const Eigen::Vector3f node2xyz[2] = {vtx2xyz.row(i_vtx0), vtx2xyz.row(i_vtx1)}; // coordinates of end points
        float w; // energy of one spring
        Eigen::Vector3f dw[2] = {Eigen::Vector3f::Zero(), Eigen::Vector3f::Zero()}; // gradient of the energy of one spring
        wdw_spring_3d( // compute energy and its gradient of a spring
            w, dw,
            node2xyz, length_ini, stiffness);
        W += w;
        // merge gradient
        for (unsigned int i_node = 0; i_node < 2; ++i_node)
        {
            int i_vtx = line2vtx(i_line, i_node);
            gradW.row(i_vtx) += dw[i_node];
        }
    }
    // adding gravitational potential energy and its gradient
    for (unsigned int i_vtx = 0; i_vtx < num_vtx; ++i_vtx)
    {
        gradW.row(i_vtx) -= mass_point * gravity;
        W -= mass_point * vtx2xyz.row(i_vtx).dot(gravity);
    }
    // adding boundary condition
    gradW = gradW.cwiseProduct(aBCFlag);
    vtx2xyz -= gradW * learning_rate;
    return W;
}


#define ACTIVATE_GAS_GEOMETRY static PRM_Name GeometryName(GAS_NAME_GEOMETRY, SIM_GEOMETRY_DATANAME); static PRM_Default GeometryNameDefault(0, SIM_GEOMETRY_DATANAME); PRMs.emplace_back(PRM_STRING, 1, &GeometryName, &GeometryNameDefault);
#define ACTIVATE_GAS_GEOMETRY_MESH static PRM_Name GeometryMeshName("geo2", "GeometryMesh"); static PRM_Default GeometryMeshNameDefault(0, "GeometryMesh"); PRMs.emplace_back(PRM_STRING, 1, &GeometryMeshName, &GeometryMeshNameDefault);
#define POINT_ATTRIBUTE_V3(NAME) GA_RWAttributeRef NAME##_attr = gdp.findGlobalAttribute(#NAME); if (!NAME##_attr.isValid()) NAME##_attr = gdp.addFloatTuple(GA_ATTRIB_POINT, #NAME, 3, GA_Defaults(0)); GA_RWHandleV3 NAME##_handle(NAME##_attr);
#define POINT_ATTRIBUTE_I(NAME) GA_RWAttributeRef NAME##_attr = gdp.findGlobalAttribute(#NAME); if (!NAME##_attr.isValid()) NAME##_attr = gdp.addIntTuple(GA_ATTRIB_POINT, #NAME, 1, GA_Defaults(0)); GA_RWHandleI NAME##_handle(NAME##_attr);
#define GLOBAL_ATTRIBUTE_F(NAME) GA_RWAttributeRef NAME##_attr = gdp.findGlobalAttribute(#NAME); if (!NAME##_attr.isValid()) NAME##_attr = gdp.addFloatTuple(GA_ATTRIB_DETAIL, #NAME, 1, GA_Defaults(0)); GA_RWHandleF NAME##_handle(NAME##_attr);

const SIM_DopDescription* GAS_Task05::getDopDescription()
{
    static std::vector<PRM_Template> PRMs;
    PRMs.clear();
    ACTIVATE_GAS_GEOMETRY
    ACTIVATE_GAS_GEOMETRY_MESH
    PRMs.emplace_back();

    static SIM_DopDescription DESC(GEN_NODE,
                                   DOP_NAME,
                                   DOP_ENGLISH,
                                   DATANAME,
                                   classname(),
                                   PRMs.data());
    DESC.setDefaultUniqueDataName(UNIQUE_DATANAME);
    setGasDescription(DESC);
    return &DESC;
}

bool GAS_Task05::solveGasSubclass(SIM_Engine& engine, SIM_Object* obj, SIM_Time time, SIM_Time timestep)
{
    SIM_GeometryCopy* G = getGeometryCopy(obj, GAS_NAME_GEOMETRY);
    SIM_GeometryAutoWriteLock lock(G);
    GU_Detail& gdp = lock.getGdp();
    POINT_ATTRIBUTE_V3(vtx2xyz_ini);
    POINT_ATTRIBUTE_I(vtx2isfree);

    SIM_GeometryCopy* Mesh = getGeometryCopy(obj, "geo2");
    SIM_GeometryAutoWriteLock lock_mesh(Mesh);
    GU_Detail& gdp_mesh = lock_mesh.getGdp();

    constexpr float learning_rate = 6.5e-4f;
    Eigen::Matrix<float, Eigen::Dynamic, 3, Eigen::RowMajor> vtx2xyz;
    Eigen::Matrix<float, Eigen::Dynamic, 3, Eigen::RowMajor> vtx2xyz_ini;
    Eigen::Matrix<int, Eigen::Dynamic, 2, Eigen::RowMajor> line2vtx;
    Eigen::MatrixX3f vtx2isfree;
    if (engine.getSimulationFrame(time) == 1)
    {
        constexpr int num_theta = 64;
        const auto [tri2vtx_ret, vtx2xyz_ini_ret] = generate_mesh_annulus3(0.3, 0.8, 32, num_theta);
        const auto line2vtx_ret = lines_of_mesh(tri2vtx_ret, static_cast<int>(vtx2xyz_ini_ret.rows()));
        vtx2xyz_ini = vtx2xyz_ini_ret;
        vtx2xyz = vtx2xyz_ini;
        line2vtx = line2vtx_ret;

        // whether the DoFs of vertices are fixed or not. Fixed: 0, Free:1
        vtx2isfree = Eigen::MatrixX3f::Ones(vtx2xyz.rows(), 3);
        for (int i = 0; i < num_theta; ++i)
        {
            vtx2isfree.row(i) = Eigen::Vector3f(0., 0., 0);
        }

        for (int i = 0; i < vtx2xyz.rows(); i++)
        {
            UT_Vector3 pos(vtx2xyz(i, 0), vtx2xyz(i, 1), vtx2xyz(i, 2));
            GA_Offset off = gdp.appendPoint();
            gdp.setPos3(off, pos);
            vtx2xyz_ini_handle.set(off, pos);

            off = gdp_mesh.appendPoint();
            gdp.setPos3(off, pos);
        }

        for (int i = 0; i < line2vtx.rows(); i++)
        {
            GA_Offset off0 = gdp.pointOffset(line2vtx(i, 0));
            GA_Offset off1 = gdp.pointOffset(line2vtx(i, 1));
            GEO_PrimPoly* poly = GU_PrimPoly::build(&gdp, 2, true, false);
            poly->setPointOffset(0, off0);
            poly->setPointOffset(1, off1);
        }

        for (int i = 0; i < num_theta; ++i)
        {
            GA_Offset off = gdp.pointOffset(i);
            vtx2isfree_handle.set(off, 1);
        }

        for (int i = 0; i < tri2vtx_ret.rows(); i++)
        {
            GA_Offset off0 = gdp_mesh.pointOffset(tri2vtx_ret(i, 0));
            GA_Offset off1 = gdp_mesh.pointOffset(tri2vtx_ret(i, 1));
            GA_Offset off2 = gdp_mesh.pointOffset(tri2vtx_ret(i, 2));
            GEO_PrimPoly* poly = GU_PrimPoly::build(&gdp_mesh, 3, false, false);
            poly->setPointOffset(0, off0);
            poly->setPointOffset(1, off2);
            poly->setPointOffset(2, off1);
        }
    }
    else
    {
        GA_Offset off;
        {
            vtx2xyz.resize(gdp.getNumPoints(), 3);
            vtx2xyz_ini.resize(gdp.getNumPoints(), 3);
            vtx2isfree = Eigen::MatrixX3f::Ones(gdp.getNumPoints(), 3);
            GA_FOR_ALL_PTOFF(&gdp, off)
            {
                UT_Vector3 pos = gdp.getPos3(off);
                GA_Index idx = gdp.pointIndex(off);
                vtx2xyz.row(idx) = Eigen::Vector3f(pos.x(), pos.y(), pos.z());
                UT_Vector3 vtx2xyz_ini_off = vtx2xyz_ini_handle.get(off);
                vtx2xyz_ini.row(idx) = Eigen::Vector3f(vtx2xyz_ini_off.x(), vtx2xyz_ini_off.y(), vtx2xyz_ini_off.z());
                int vtx2isfree_off = vtx2isfree_handle.get(off);
                if (vtx2isfree_off == 1)
                    vtx2isfree.row(idx) = Eigen::Vector3f(0., 0., 0);
            }
        }

        {
            const GEO_Primitive* prim;
            line2vtx.resize(gdp.getNumPrimitives(), 2);
            int iter = 0;
            GA_FOR_ALL_PRIMITIVES(&gdp, prim)
            {
                GA_Offset off0 = prim->getPointOffset(0);
                GA_Offset off1 = prim->getPointOffset((1));
                GA_Index idx0 = gdp.pointIndex(off0);
                GA_Index idx1 = gdp.pointIndex(off1);
                line2vtx.row(iter++) = Eigen::Vector2i(idx0, idx1);
            }
        }
    }

    for (int itr = 0; itr < 400; ++itr)
    {
        float W = gradient_descent_energy_minimization(
            vtx2xyz, vtx2xyz_ini, line2vtx, 60.f, 1.f, {0., -0.1, 0}, vtx2isfree, learning_rate);
    }

    {
        GA_Offset off;
        GA_FOR_ALL_PTOFF(&gdp, off)
        {
            GA_Index idx = gdp.pointIndex(off);
            UT_Vector3 pos(vtx2xyz(idx, 0), vtx2xyz(idx, 1), vtx2xyz(idx, 2));
            gdp.setPos3(off, pos);
        }
    }

    {
        GA_Offset off;
        GA_FOR_ALL_PTOFF(&gdp_mesh, off)
        {
            GA_Index idx = gdp_mesh.pointIndex(off);
            UT_Vector3 pos(vtx2xyz(idx, 0), vtx2xyz(idx, 1), vtx2xyz(idx, 2));
            gdp_mesh.setPos3(off, pos);
        }
    }

    return true;
}
