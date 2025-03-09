#include "Task10.h"

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

/**
 * inertia tensor of 3D solid triangle mesh
 * @param tri2vtx connectivity of a triangle mesh
 * @param vtx2xyz coordinates
 * @return
 */
Eigen::Matrix3f inertia_tensor_solid_3d_triangle_mesh(
    const Eigen::Matrix<int, Eigen::Dynamic, 3, Eigen::RowMajor>& tri2vtx,
    const Eigen::Matrix<float, Eigen::Dynamic, 3, Eigen::RowMajor>& vtx2xyz)
{
    Eigen::Matrix3f res = Eigen::Matrix3f::Zero();
    for (int i_tri = 0; i_tri < tri2vtx.rows(); ++i_tri)
    {
        const Eigen::Vector3f p0 = vtx2xyz.row(tri2vtx(i_tri, 0));
        const Eigen::Vector3f p1 = vtx2xyz.row(tri2vtx(i_tri, 1));
        const Eigen::Vector3f p2 = vtx2xyz.row(tri2vtx(i_tri, 2));
        Eigen::Vector3f pa = p0 + p1 + p2;
        const Eigen::Matrix3f m0 = p0 * p0.transpose();
        const Eigen::Matrix3f m1 = p1 * p1.transpose();
        const Eigen::Matrix3f m2 = p2 * p2.transpose();
        const Eigen::Matrix3f ma = pa * pa.transpose();
        const Eigen::Matrix3f I0 = m0 + m1 + m2 + ma;
        float tr0 = I0.trace();
        Eigen::Matrix3f I = tr0 * Eigen::Matrix3f::Identity() - I0;
        float v = p0.dot(p1.cross(p2));
        res += (v / 120.f) * I;
    }
    return res;
}

inline Eigen::Matrix3f skew(const Eigen::Vector3f& v)
{
    Eigen::Matrix3f S;
    S << 0.0f, -v.z(), v.y(),
        v.z(), 0.0f, -v.x(),
        -v.y(), v.x(), 0.0f;
    return S;
}

#define ACTIVATE_GAS_GEOMETRY static PRM_Name GeometryName(GAS_NAME_GEOMETRY, SIM_GEOMETRY_DATANAME); static PRM_Default GeometryNameDefault(0, SIM_GEOMETRY_DATANAME); PRMs.emplace_back(PRM_STRING, 1, &GeometryName, &GeometryNameDefault);
#define POINT_ATTRIBUTE_V3(NAME) GA_RWAttributeRef NAME##_attr = gdp.findGlobalAttribute(#NAME); if (!NAME##_attr.isValid()) NAME##_attr = gdp.addFloatTuple(GA_ATTRIB_POINT, #NAME, 3, GA_Defaults(0)); GA_RWHandleV3 NAME##_handle(NAME##_attr);

const SIM_DopDescription* GAS_Task10::getDopDescription()
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

void GAS_Task10::initializeSubclass()
{
    GAS_SubSolver::initializeSubclass();
    this->Omega = Eigen::Vector3f(0.f, 0.05f, 1.f); // initial angular velocity (\dot{R} = R * Skew(\Omega))
    this->rotation = Eigen::Matrix3f::Identity(); // rotation to optimize
}

void GAS_Task10::makeEqualSubclass(const SIM_Data* source)
{
    GAS_SubSolver::makeEqualSubclass(source);
    this->Omega = ((GAS_Task10*)source)->Omega;
    this->rotation = ((GAS_Task10*)source)->rotation;
}

bool GAS_Task10::solveGasSubclass(SIM_Engine& engine, SIM_Object* obj, SIM_Time time, SIM_Time timestep)
{
    SIM_GeometryCopy* G = getGeometryCopy(obj, GAS_NAME_GEOMETRY);
    SIM_GeometryAutoWriteLock lock(G);
    GU_Detail& gdp = lock.getGdp();
    POINT_ATTRIBUTE_V3(vtx2xyz_ini);

    Eigen::Matrix<float, Eigen::Dynamic, 3, Eigen::RowMajor> vtx2xyz;
    Eigen::Matrix<float, Eigen::Dynamic, 3, Eigen::RowMajor> vtx2xyz_ini;
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

    const Eigen::Matrix3f inertia = inertia_tensor_solid_3d_triangle_mesh(tri2vtx, vtx2xyz_ini);


    constexpr float dt = 0.001; // time step
    constexpr unsigned int i_vtx_trajectory = 831;

    for (int itr = 0; itr < 100; ++itr)
    {
        // sub-stepping
        // Write some code below to simulate rotation of the rigid body
        // Use the **forward Euler method** to update the rotation matrix and the angular velocity
        // Note that `rotation` should stay a rotational matrix after the update
        rotation = rotation * (Eigen::Matrix3f::Identity() + dt * skew(Omega));
        Omega = Omega + dt * inertia.inverse() * ((inertia * Omega).cross(Omega));
        // Do not change anything else except for the two lines above.
    }
    vtx2xyz = (rotation * vtx2xyz_ini.transpose()).transpose(); // the rotated mesh's vertices


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
