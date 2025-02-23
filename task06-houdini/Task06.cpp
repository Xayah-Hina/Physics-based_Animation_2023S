#include "Task06.h"

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

    auto vertex_to_vertex(
        const Eigen::MatrixXi& elem2vtx,
        size_t num_vtx)
    {
        const auto [vtx2idx, idx2elem] = vertex_to_elem(elem2vtx, num_vtx);
        std::vector<unsigned int> vtx2jdx, jdx2vtx;
        vtx2jdx.push_back(0);
        for (unsigned int i_vtx = 0; i_vtx < num_vtx; ++i_vtx)
        {
            std::set<int> set_connected_points;
            for (unsigned int idx = vtx2idx[i_vtx]; idx < vtx2idx[i_vtx + 1]; ++idx)
            {
                const unsigned int i_elem = idx2elem[idx];
                for (int i_node = 0; i_node < elem2vtx.cols(); ++i_node)
                {
                    const int j_vtx = elem2vtx(i_elem, i_node);
                    set_connected_points.insert(j_vtx);
                }
            }
            vtx2jdx.push_back(vtx2jdx[vtx2jdx.size() - 1] + set_connected_points.size());
            jdx2vtx.insert(jdx2vtx.end(), set_connected_points.begin(), set_connected_points.end());
        }
        assert(vtx2jdx.size() == num_vtx+1);
        assert(vtx2jdx[num_vtx] == jdx2vtx.size());
        return std::make_pair(vtx2jdx, jdx2vtx);
    }

    template <int N>
    class BlockSparseMatrix
    {
        using Vector = Eigen::Matrix<double, Eigen::Dynamic, N>;
        using BlockVector = Eigen::Vector<double, N>;
        using BlockMatrix = Eigen::Matrix<double, N, N>;

    public:
        void initialize(const Eigen::MatrixXi& elem2vtx, unsigned int num_vtx)
        {
            auto [vtx2idx, idx2vtx] = pba::vertex_to_vertex(elem2vtx, num_vtx);
            this->row2idx = vtx2idx;
            this->idx2col = idx2vtx;
            idx2block.resize(this->idx2col.size());
        }

        void setZero()
        {
            for (auto& block : idx2block)
            {
                block.setZero();
            }
        }

        auto& coeff(unsigned int i_row, unsigned int i_col)
        {
            auto itr0 = idx2col.begin() + row2idx[i_row];
            auto itr1 = idx2col.begin() + row2idx[i_row + 1];
            auto itr2 = std::find(itr0, itr1, i_col);
            assert(itr2 != itr1);
            unsigned int idx0 = std::distance(idx2col.begin(), itr2);
            return idx2block[idx0];
        }

        void set_is_free(
            const Vector& vtx2isfree)
        {
            unsigned int num_row = row2idx.size() - 1;
            for (unsigned int irow = 0; irow < num_row; ++irow)
            {
                for (unsigned int idx0 = row2idx[irow]; idx0 < row2idx[irow + 1]; ++idx0)
                {
                    unsigned int icol = idx2col[idx0];
                    BlockMatrix dia_row = vtx2isfree.row(irow).asDiagonal();
                    BlockMatrix dia_col = vtx2isfree.row(icol).asDiagonal();
                    idx2block[idx0] = dia_col * idx2block[idx0] * dia_row;
                    if (icol == irow)
                    {
                        for (int i = 0; i < N; ++i)
                        {
                            idx2block[idx0](i, i) += 1.0 - vtx2isfree(irow, i);
                        }
                    }
                }
            }
        }

        Vector solve_conjugate_gradient(
            Vector& r) const
        {
            Vector x = Eigen::MatrixX3d::Zero(r.rows(), r.cols());
            Vector p = r;
            Vector Ap = p;
            const double r_squared_norm_ini = r.squaredNorm();
            double r_squared_norm_pre = r_squared_norm_ini;
            for (unsigned int itr = 0; itr < 10; ++itr)
            {
                this->multiply_vector(Ap, p);
                double alpha = r_squared_norm_pre / (p.cwiseProduct(Ap)).sum();
                x += alpha * p;
                r -= alpha * Ap;
                const double r_squared_norm_pos = r.squaredNorm();
                if (r_squared_norm_pos < r_squared_norm_ini * 1.0e-4) { return x; }
                double beta = r_squared_norm_pos / r_squared_norm_pre;
                r_squared_norm_pre = r_squared_norm_pos;
                p = r + beta * p;
            }
            return x;
        }

    private:
        void multiply_vector(Vector& y,
                             const Vector& x) const
        {
            unsigned int num_row = row2idx.size() - 1;
            y.setZero();
            for (unsigned int irow = 0; irow < num_row; ++irow)
            {
                for (unsigned int idx0 = row2idx[irow]; idx0 < row2idx[irow + 1]; ++idx0)
                {
                    unsigned int icol = idx2col[idx0];
                    BlockVector x0 = x.row(icol);
                    y.row(irow) += idx2block[idx0] * x0;
                }
            }
        }

    public:
        std::vector<unsigned int> row2idx;
        std::vector<unsigned int> idx2col;
        std::vector<BlockMatrix> idx2block;
    };

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
}

/**
 * compute the elastic potential energy, its gradient and its hessian of a 3D spring.
 * @param w elastic potential energy
 * @param dw gradient of energy
 * @param ddw hessian of energy
 * @param node2xyz xyz coordinates of two end points
 * @param length_ini initial length
 * @param stiffness Hock's spring coefficient
 */
void WdWddW_Spring3(
    double& w,
    Eigen::Vector3d dw[2],
    Eigen::Matrix3d ddw[2][2],
    const Eigen::Vector3d node2xyz[2],
    const double length_ini,
    const double stiffness)
{
    constexpr int num_node = 2; // number of two end points
    const double length = (node2xyz[0] - node2xyz[1]).norm(); // distance between p0 and p1
    const double C = length - length_ini; // the length differences.
    w = 0.5f * stiffness * C * C; // Hooke's law. energy is square of length difference W=1/2*k*C*C
    //
    const Eigen::Vector3d u01 = (node2xyz[1] - node2xyz[0]).normalized();
    const Eigen::Vector3d dC[num_node] = {-u01, u01};
    for (int ino = 0; ino < num_node; ++ino)
    {
        dw[ino] = stiffness * dC[ino] * C; // dW = k*dC*C
    }
    // write code to correctly compute hessian of a spring.
    // ddw[i_node][j_node] stands for derivative of dw[i_node] w.r.t the end-point's position node2xyz[j_node]
    // the current hessian computed by the code below is not very accurate, so the simulation is unstable.
    // const Eigen::Matrix3d n = stiffness * u01 * u01.transpose();
    // ddw[0][0] = n;
    // ddw[1][1] = n;
    // ddw[0][1] = -n;
    // ddw[1][0] = -n;

    // 计算单位矩阵
    const Eigen::Matrix3d I = Eigen::Matrix3d::Identity();
    // 计算矩阵 A
    Eigen::Matrix3d A = u01 * u01.transpose() + (C / length) * (I - u01 * u01.transpose());
    // 填充 Hessian (ddw[i][j] 表示 dw[i] 关于 node2xyz[j] 的二阶导数)
    ddw[0][0] = stiffness * A;
    ddw[0][1] = -stiffness * A;
    ddw[1][0] = -stiffness * A;
    ddw[1][1] = stiffness * A;
}

float step_time_mass_spring_system_with_variational_integration(
    Eigen::Matrix<float, Eigen::Dynamic, 3, Eigen::RowMajor>& vtx2xyz,
    Eigen::Matrix<float, Eigen::Dynamic, 3, Eigen::RowMajor>& vtx2velocity,
    const Eigen::Matrix<float, Eigen::Dynamic, 3, Eigen::RowMajor>& vtx2xyz_ini,
    const Eigen::Matrix<int, Eigen::Dynamic, 2, Eigen::RowMajor>& line2vtx,
    float stiffness,
    float mass_point,
    const Eigen::Vector3f& gravity,
    const Eigen::MatrixX3d& vtx2isfree,
    float dt,
    pba::BlockSparseMatrix<3>& sparse)
{
    // simulation
    const unsigned int num_vtx = vtx2xyz.rows(); // number of vertices
    double W = 0.0; // energy of the system
    Eigen::MatrixX3d gradW = Eigen::MatrixX3d::Zero(num_vtx, 3); // gradient of the energy
    sparse.setZero();
    // step position using velocity
    vtx2xyz += vtx2velocity * dt;
    for (int i_line = 0; i_line < line2vtx.rows(); ++i_line)
    {
        // loop over springs
        const int i_vtx0 = line2vtx(i_line, 0); // index of vertex 0
        const int i_vtx1 = line2vtx(i_line, 1); // index of vertex 1
        const double length_ini = (vtx2xyz_ini.row(i_vtx0) - vtx2xyz_ini.row(i_vtx1)).norm(); // initial length
        const Eigen::Vector3d node2xyz[2] = {
            // coordinates of end points
            vtx2xyz.row(i_vtx0).cast<double>(),
            vtx2xyz.row(i_vtx1).cast<double>()
        };
        double w; // energy of one spring
        Eigen::Vector3d dw[2]; // gradient of the energy of one spring
        Eigen::Matrix3d ddw[2][2]; // hessian of energy of one spring
        WdWddW_Spring3(
            w, dw, ddw,
            node2xyz, length_ini, stiffness);
        // merge energy
        W += w;
        // merge gradient
        for (unsigned int i_node = 0; i_node < 2; ++i_node)
        {
            const int i_vtx = line2vtx(i_line, i_node);
            gradW.row(i_vtx) += dw[i_node];
        }
        // merge hessian
        for (unsigned int i_node = 0; i_node < 2; ++i_node)
        {
            for (unsigned int j_node = 0; j_node < 2; ++j_node)
            {
                const int i_vtx = line2vtx(i_line, i_node);
                const int j_vtx = line2vtx(i_line, j_node);
                sparse.coeff(i_vtx, j_vtx) += ddw[i_node][j_node];
            }
        }
    }
    // adding the dynamic effect
    for (unsigned int i_vtx = 0; i_vtx < num_vtx; ++i_vtx)
    {
        sparse.coeff(i_vtx, i_vtx) += Eigen::Matrix3d::Identity() * (mass_point / (dt * dt));
    }
    // adding gravitational potential energy and its gradient
    for (unsigned int i_vtx = 0; i_vtx < num_vtx; ++i_vtx)
    {
        gradW.row(i_vtx) -= (mass_point * gravity).cast<double>();
        W -= mass_point * vtx2xyz.row(i_vtx).dot(gravity);
    }
    // set free/fix
    gradW = gradW.cwiseProduct(vtx2isfree);
    sparse.set_is_free(vtx2isfree);
    // solve
    Eigen::MatrixX3d x = sparse.solve_conjugate_gradient(gradW);
    // step position and velocity
    vtx2velocity += -x.cast<float>() / dt;
    vtx2xyz -= x.cast<float>();
    return static_cast<float>(W);
}


#define ACTIVATE_GAS_GEOMETRY static PRM_Name GeometryName(GAS_NAME_GEOMETRY, SIM_GEOMETRY_DATANAME); static PRM_Default GeometryNameDefault(0, SIM_GEOMETRY_DATANAME); PRMs.emplace_back(PRM_STRING, 1, &GeometryName, &GeometryNameDefault);
#define ACTIVATE_GAS_GEOMETRY_MESH static PRM_Name GeometryMeshName("geo2", "GeometryMesh"); static PRM_Default GeometryMeshNameDefault(0, "GeometryMesh"); PRMs.emplace_back(PRM_STRING, 1, &GeometryMeshName, &GeometryMeshNameDefault);
#define POINT_ATTRIBUTE_V3(NAME) GA_RWAttributeRef NAME##_attr = gdp.findGlobalAttribute(#NAME); if (!NAME##_attr.isValid()) NAME##_attr = gdp.addFloatTuple(GA_ATTRIB_POINT, #NAME, 3, GA_Defaults(0)); GA_RWHandleV3 NAME##_handle(NAME##_attr);
#define POINT_ATTRIBUTE_I(NAME) GA_RWAttributeRef NAME##_attr = gdp.findGlobalAttribute(#NAME); if (!NAME##_attr.isValid()) NAME##_attr = gdp.addIntTuple(GA_ATTRIB_POINT, #NAME, 1, GA_Defaults(0)); GA_RWHandleI NAME##_handle(NAME##_attr);
#define GLOBAL_ATTRIBUTE_F(NAME) GA_RWAttributeRef NAME##_attr = gdp.findGlobalAttribute(#NAME); if (!NAME##_attr.isValid()) NAME##_attr = gdp.addFloatTuple(GA_ATTRIB_DETAIL, #NAME, 1, GA_Defaults(0)); GA_RWHandleF NAME##_handle(NAME##_attr);

const SIM_DopDescription* GAS_Task06::getDopDescription()
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

bool GAS_Task06::solveGasSubclass(SIM_Engine& engine, SIM_Object* obj, SIM_Time time, SIM_Time timestep)
{
    SIM_GeometryCopy* G = getGeometryCopy(obj, GAS_NAME_GEOMETRY);
    SIM_GeometryAutoWriteLock lock(G);
    GU_Detail& gdp = lock.getGdp();
    POINT_ATTRIBUTE_V3(vtx2xyz_ini);
    POINT_ATTRIBUTE_V3(vtx2velocity);
    POINT_ATTRIBUTE_I(vtx2isfree);

    SIM_GeometryCopy* Mesh = getGeometryCopy(obj, "geo2");
    SIM_GeometryAutoWriteLock lock_mesh(Mesh);
    GU_Detail& gdp_mesh = lock_mesh.getGdp();

    Eigen::Matrix<float, Eigen::Dynamic, 3, Eigen::RowMajor> vtx2xyz;
    Eigen::Matrix<float, Eigen::Dynamic, 3, Eigen::RowMajor> vtx2xyz_ini;
    Eigen::Matrix<float, Eigen::Dynamic, 3, Eigen::RowMajor> vtx2velocity;
    Eigen::Matrix<int, Eigen::Dynamic, 2, Eigen::RowMajor> line2vtx;
    Eigen::Matrix<int, Eigen::Dynamic, 3, Eigen::RowMajor> tri2vtx;
    Eigen::MatrixX3d vtx2isfree;

    if (engine.getSimulationFrame(time) == 1)
    {
        constexpr int num_theta = 64;
        const auto [tri2vtx_ret, vtx2xyz_ini_ret] = pba::generate_mesh_annulus3(0.3, 0.8, 32, num_theta);
        const auto line2vtx_ret = pba::lines_of_mesh(tri2vtx_ret, static_cast<int>(vtx2xyz_ini_ret.rows()));
        vtx2xyz_ini = vtx2xyz_ini_ret;
        vtx2xyz = vtx2xyz_ini;
        line2vtx = line2vtx_ret;
        tri2vtx = tri2vtx_ret;

        vtx2velocity.resize(vtx2xyz.rows(), 3);
        vtx2velocity.setZero();

        // whether the DoFs of vertices are fixed or not. Fixed: 0, Free:1
        vtx2isfree = Eigen::MatrixX3d::Ones(vtx2xyz.rows(), 3);
        for (int i = 0; i < num_theta; ++i)
        {
            vtx2isfree.row(i) = Eigen::Vector3d(0., 0., 0);
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
            vtx2isfree = Eigen::MatrixX3d::Ones(vtx2xyz.rows(), 3);
            vtx2velocity.resize(gdp.getNumPoints(), 3);
            GA_FOR_ALL_PTOFF(&gdp, off)
            {
                UT_Vector3 pos = gdp.getPos3(off);
                GA_Index idx = gdp.pointIndex(off);
                vtx2xyz.row(idx) = Eigen::Vector3f(pos.x(), pos.y(), pos.z());
                UT_Vector3 vtx2xyz_ini_off = vtx2xyz_ini_handle.get(off);
                vtx2xyz_ini.row(idx) = Eigen::Vector3f(vtx2xyz_ini_off.x(), vtx2xyz_ini_off.y(), vtx2xyz_ini_off.z());
                int vtx2isfree_off = vtx2isfree_handle.get(off);
                if (vtx2isfree_off == 1)
                    vtx2isfree.row(idx) = Eigen::Vector3d(0., 0., 0);
                UT_Vector3 vel = vtx2velocity_handle.get(off);
                vtx2velocity.row(idx) = Eigen::Vector3f(vel.x(), vel.y(), vel.z());
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

        {
            const GEO_Primitive* prim;
            tri2vtx.resize(gdp_mesh.getNumPrimitives(), 3);
            int iter = 0;
            GA_FOR_ALL_PRIMITIVES(&gdp_mesh, prim)
            {
                GA_Offset off0 = prim->getPointOffset(0);
                GA_Offset off1 = prim->getPointOffset(1);
                GA_Offset off2 = prim->getPointOffset(2);
                GA_Index idx0 = gdp_mesh.pointIndex(off0);
                GA_Index idx1 = gdp_mesh.pointIndex(off1);
                GA_Index idx2 = gdp_mesh.pointIndex(off2);
                tri2vtx.row(iter++) = Eigen::Vector3i(idx0, idx2, idx1);
            }
        }
    }

    const float dt = 0.13f;
    pba::BlockSparseMatrix<3> sparse_matrix;
    sparse_matrix.initialize(tri2vtx, vtx2xyz.rows());
    float W = step_time_mass_spring_system_with_variational_integration(
        vtx2xyz, vtx2velocity, vtx2xyz_ini, line2vtx, 60.f, 1.f, {0., -0.1, 0}, vtx2isfree, dt,
        sparse_matrix);

    std::cout << "W: " << W << std::endl;

    {
        GA_Offset off;
        GA_FOR_ALL_PTOFF(&gdp, off)
        {
            GA_Index idx = gdp.pointIndex(off);
            UT_Vector3 pos(vtx2xyz(idx, 0), vtx2xyz(idx, 1), vtx2xyz(idx, 2));
            UT_Vector3 vel(vtx2velocity(idx, 0), vtx2velocity(idx, 1), vtx2velocity(idx, 2));
            gdp.setPos3(off, pos);
            vtx2velocity_handle.set(off, vel);
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
