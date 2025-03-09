#include "Task08.h"

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

#ifndef M_PI
  #define M_PI 3.14159265358979323846264338327950288
#endif

namespace pba
{
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
}

/**
 * compute elastic potential energy, energy's gradient, and energy's hessian
 * @param w energy
 * @param dw gradient
 * @param ddw hessian
 * @param node2xyz node's xyz-coordinates
 * @param length_ini initial length
 * @param stiffness spring's stiffness parameter
 */
void wdwddw_spring(
    double& w,
    Eigen::Vector3d dw[2],
    Eigen::Matrix3d ddw[2][2],
    const Eigen::Vector3d node2xyz[2],
    double length_ini,
    double stiffness)
{
    double length = (node2xyz[1] - node2xyz[0]).norm();
    double C = length - length_ini;
    w = 0.5 * stiffness * C * C;
    Eigen::Vector3d dC = (node2xyz[0] - node2xyz[1]).normalized();
    Eigen::Matrix3d ddC = (1. / length) * (Eigen::Matrix3d::Identity() - dC * dC.transpose());
    dw[0] = stiffness * C * dC;
    dw[1] = -dw[0];
    ddw[0][0] = stiffness * dC * dC.transpose() + stiffness * C * ddC;
    ddw[1][1] = ddw[0][0];
    ddw[1][0] = -ddw[0][0];
    ddw[0][1] = -ddw[0][0];
}

/**
 * compute volume of a tetrahedron connecting vertices of the triangle and the origin
 * @param w volume
 * @param dw differentiation of volume w.r.t the nodes' xyz coordinates
 * @param node2xyz nodes' xyz coordinates
 */
void wdw_volume_tri_origin(
    double& w,
    Eigen::Vector3d dw[3],
    const Eigen::Vector3d node2xyz[3])
{
    w = 1. / 6. * node2xyz[0].dot(node2xyz[1].cross(node2xyz[2]));
    // ------------------------------
    // Write some code below to compute differentiation. Keep it simple and understandable
    dw[0] = (1.0 / 6.0) * (node2xyz[1].cross(node2xyz[2]));
    dw[1] = (1.0 / 6.0) * (node2xyz[2].cross(node2xyz[0]));
    dw[2] = (1.0 / 6.0) * (node2xyz[0].cross(node2xyz[1]));
}

void inflate(
    Eigen::Matrix<float, Eigen::Dynamic, 3, Eigen::RowMajor>& vtx2xyz,
    double& lambda,
    double volume_trg,
    const Eigen::Matrix<int, Eigen::Dynamic, 3, Eigen::RowMajor>& tri2vtx,
    const Eigen::Matrix<int, Eigen::Dynamic, 2, Eigen::RowMajor>& line2vtx,
    const Eigen::Matrix<float, Eigen::Dynamic, 3, Eigen::RowMajor>& vtx2xyz_ini)
{
    double stiffness = 1.0;
    unsigned int num_vtx = vtx2xyz.rows();
    Eigen::MatrixXd ddW(num_vtx * 3 + 1, num_vtx * 3 + 1);
    Eigen::VectorXd dW(num_vtx * 3 + 1);
    ddW.setZero();
    dW.setZero();
    // setting mass-spring system as a regularizer
    double elastic_energy = 0.0;
    for (unsigned int i_line = 0; i_line < line2vtx.rows(); ++i_line)
    {
        const int node2vtx[2] = {
            // index of two end points
            line2vtx(i_line, 0),
            line2vtx(i_line, 1)
        };
        float length_ini = (vtx2xyz_ini.row(node2vtx[0]) - vtx2xyz_ini.row(node2vtx[1])).norm();
        const Eigen::Vector3d node2xyz[3] = {
            // xyz-coordinates of two nodes
            vtx2xyz.row(node2vtx[0]).cast<double>(),
            vtx2xyz.row(node2vtx[1]).cast<double>()
        };
        double w = 0;
        Eigen::Vector3d dw[2];
        Eigen::Matrix3d ddw[2][2];
        wdwddw_spring(w, dw, ddw,
                      node2xyz, length_ini, stiffness);
        elastic_energy += w;
        for (unsigned int inode = 0; inode < 2; ++inode)
        {
            for (unsigned int idim = 0; idim < 3; ++idim)
            {
                for (unsigned int jnode = 0; jnode < 2; ++jnode)
                {
                    for (unsigned int jdim = 0; jdim < 3; ++jdim)
                    {
                        ddW(node2vtx[inode] * 3 + idim, node2vtx[jnode] * 3 + jdim) += ddw[inode][jnode](idim, jdim);
                    }
                }
                dW(node2vtx[inode] * 3 + idim) += dw[inode](idim);
            }
        }
    }
    // setting constraint
    double volume = 0.0;
    for (unsigned int i_tri = 0; i_tri < tri2vtx.rows(); ++i_tri)
    {
        const int node2vtx[3] = {
            tri2vtx(i_tri, 0),
            tri2vtx(i_tri, 1),
            tri2vtx(i_tri, 2)
        };
        const Eigen::Vector3d node2xyz[3] = {
            vtx2xyz.row(node2vtx[0]).cast<double>(),
            vtx2xyz.row(node2vtx[1]).cast<double>(),
            vtx2xyz.row(node2vtx[2]).cast<double>()
        };
        double w = 0.0;
        Eigen::Vector3d dw[3];
        wdw_volume_tri_origin(
            w, dw,
            node2xyz);
        volume += w;
        // -------------------------
        // write some code below to set values in the linear system to set constraint to specify volume
        // write some code including 'w'
        for (unsigned int inode = 0; inode < 3; ++inode)
        {
            for (unsigned int idim = 0; idim < 3; ++idim)
            {
                // Add derivative of the constraint to dW
                dW(node2vtx[inode] * 3 + idim) += lambda * dw[inode](idim);

                // Coupling in the Hessian for x–λ and λ–x
                ddW(node2vtx[inode] * 3 + idim, num_vtx * 3) += dw[inode](idim);
                ddW(num_vtx * 3, node2vtx[inode] * 3 + idim) += dw[inode](idim);
            }
        }
    }
    // Do not forget to write one line of code here
    // -------------------------------------------------
    dW(num_vtx * 3) += (volume - volume_trg);
    // Do not change below
    // damping for stable convergence
    for (unsigned int i = 0; i < ddW.rows(); ++i)
    {
        ddW(i, i) += 0.1;
    }
    std::cout << "   elastic_energy (write down in the README.md): " << elastic_energy << std::endl;
    std::cout << "   current volume: " << volume << "  target volume: " << volume_trg << std::endl;
    std::cout << "   residual (make sure this number get smaller in each iteration): " << dW.squaredNorm() << std::endl;
    // solve linear system
    Eigen::VectorXd upd = ddW.inverse() * dW;
    // update solution
    for (unsigned int i_vtx = 0; i_vtx < num_vtx; ++i_vtx)
    {
        for (unsigned int idim = 0; idim < 3; ++idim)
        {
            vtx2xyz(i_vtx, idim) -= static_cast<float>(upd(i_vtx * 3 + idim));
        }
    }
    lambda -= upd(num_vtx * 3);
    std::cout << "   lambda (make sure this number converges): " << lambda << std::endl;
}

#define ACTIVATE_GAS_GEOMETRY static PRM_Name GeometryName(GAS_NAME_GEOMETRY, SIM_GEOMETRY_DATANAME); static PRM_Default GeometryNameDefault(0, SIM_GEOMETRY_DATANAME); PRMs.emplace_back(PRM_STRING, 1, &GeometryName, &GeometryNameDefault);
#define GLOBAL_ATTRIBUTE_F(NAME) GA_RWAttributeRef NAME##_attr = gdp.findGlobalAttribute(#NAME); if (!NAME##_attr.isValid()) NAME##_attr = gdp.addFloatTuple(GA_ATTRIB_DETAIL, #NAME, 1, GA_Defaults(0)); GA_RWHandleF NAME##_handle(NAME##_attr);
#define POINT_ATTRIBUTE_V3(NAME) GA_RWAttributeRef NAME##_attr = gdp.findGlobalAttribute(#NAME); if (!NAME##_attr.isValid()) NAME##_attr = gdp.addFloatTuple(GA_ATTRIB_POINT, #NAME, 3, GA_Defaults(0)); GA_RWHandleV3 NAME##_handle(NAME##_attr);

const SIM_DopDescription* GAS_Task08::getDopDescription()
{
    static std::vector<PRM_Template> PRMs;
    PRMs.clear();
    ACTIVATE_GAS_GEOMETRY
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

bool GAS_Task08::solveGasSubclass(SIM_Engine& engine, SIM_Object* obj, SIM_Time time, SIM_Time timestep)
{
    SIM_GeometryCopy* G = getGeometryCopy(obj, GAS_NAME_GEOMETRY);
    SIM_GeometryAutoWriteLock lock(G);
    GU_Detail& gdp = lock.getGdp();
    POINT_ATTRIBUTE_V3(vtx2xyz_ini);
    GLOBAL_ATTRIBUTE_F(lambda);
    GLOBAL_ATTRIBUTE_F(volume_trg);

    Eigen::Matrix<float, Eigen::Dynamic, 3, Eigen::RowMajor> vtx2xyz;
    Eigen::Matrix<float, Eigen::Dynamic, 3, Eigen::RowMajor> vtx2xyz_ini;
    Eigen::Matrix<int, Eigen::Dynamic, 2, Eigen::RowMajor> line2vtx;
    Eigen::Matrix<int, Eigen::Dynamic, 3, Eigen::RowMajor> tri2vtx;


    vtx2xyz.resize(gdp.getNumPoints(), 3);
    vtx2xyz_ini.resize(gdp.getNumPoints(), 3);
    tri2vtx.resize(gdp.getNumPrimitives(), 3);

    GA_Offset off;
    {
        GA_FOR_ALL_PTOFF(&gdp, off)
        {
            UT_Vector3 pos = gdp.getPos3(off);
            GA_Index idx = gdp.pointIndex(off);
            vtx2xyz.row(idx) = Eigen::Vector3f(pos.x(), pos.y(), pos.z());
            UT_Vector3 vtx2xyz_ini_off = vtx2xyz_ini_handle.get(off);
            vtx2xyz_ini.row(idx) = Eigen::Vector3f(vtx2xyz_ini_off.x(), vtx2xyz_ini_off.y(), vtx2xyz_ini_off.z());
        }
    }
    {
        const GEO_Primitive* prim;
        int iter = 0;
        GA_FOR_ALL_PRIMITIVES(&gdp, prim)
        {
            GA_Offset off0 = prim->getPointOffset(0);
            GA_Offset off1 = prim->getPointOffset(1);
            GA_Offset off2 = prim->getPointOffset(2);
            GA_Index idx0 = gdp.pointIndex(off0);
            GA_Index idx1 = gdp.pointIndex(off1);
            GA_Index idx2 = gdp.pointIndex(off2);
            tri2vtx.row(iter++) = Eigen::Vector3i(idx0, idx2, idx1);
        }
    }
    line2vtx = pba::lines_of_mesh(tri2vtx, static_cast<int>(vtx2xyz_ini.rows()));

    if (engine.getSimulationFrame(time) == 1)
    {
        double lambda = 0.0;
        double volume_ini = 0.0;
        for (unsigned int i_tri = 0; i_tri < tri2vtx.rows(); ++i_tri)
        {
            const Eigen::Vector3d node2xyz[3] = {
                vtx2xyz.row(tri2vtx(i_tri, 0)).cast<double>(),
                vtx2xyz.row(tri2vtx(i_tri, 1)).cast<double>(),
                vtx2xyz.row(tri2vtx(i_tri, 2)).cast<double>()
            };
            double w = 0.0;
            Eigen::Vector3d dw[3];
            wdw_volume_tri_origin(
                w, dw,
                node2xyz);
            volume_ini += w;
        }
        const double volume_trg = volume_ini * 2.0;

        lambda_handle.set(0, lambda);
        volume_trg_handle.set(0, volume_trg);
    }

    double lambda = lambda_handle.get(0);
    const double volume_trg = volume_trg_handle.get(0);

    inflate(vtx2xyz, lambda, volume_trg, tri2vtx, line2vtx, vtx2xyz_ini);

    lambda_handle.set(0, lambda);
    volume_trg_handle.set(0, volume_trg);

    {
        GA_FOR_ALL_PTOFF(&gdp, off)
        {
            GA_Index idx = gdp.pointIndex(off);
            UT_Vector3 pos(vtx2xyz(idx, 0), vtx2xyz(idx, 1), vtx2xyz(idx, 2));
            gdp.setPos3(off, pos);
        }
    }

    return true;
}
