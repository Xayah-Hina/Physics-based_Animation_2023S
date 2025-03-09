#include "Task09.h"

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
 *
 * @param [out] rotation rotation matrix
 * @param [out] translation translation vector
 * @param [out] vtx2xyz the transformed coordinates of the mesh's vertices
 * @param [in] i_vtx_fix the index of the fixed vertex
 * @param [in] vtx2xyz_ini the initial coordinates of the mesh's vertices
 */
void step(
    Eigen::Matrix3f& rotation,
    Eigen::Vector3f& translation,
    Eigen::Matrix<float, Eigen::Dynamic, 3, Eigen::RowMajor>& vtx2xyz,
    unsigned int i_vtx_fix,
    const Eigen::Matrix<float, Eigen::Dynamic, 3, Eigen::RowMajor>& vtx2xyz_ini)
{
    Eigen::Vector3f gravity(0.f, -10.f, 0.f);
    constexpr float penalty = 1.0e+6;
    constexpr float learning_rate = 1.0e-6;
    //
    vtx2xyz = ((rotation * vtx2xyz_ini.transpose()).colwise() + translation).transpose();
    // energy = 0.5*stiffness*||[R]*{X_fix}+{t}-{X_fix}||^2 - sum<{gravity}^T([R]*{X}+{t})>
    float E = 0.5f * penalty * (vtx2xyz.row(i_vtx_fix) - vtx2xyz_ini.row(i_vtx_fix)).squaredNorm() - (vtx2xyz * gravity).sum();
    std::cout << "energy: " << E << std::endl;
    Eigen::Vector3f dEdo = Eigen::Vector3f::Zero();
    Eigen::Vector3f dEdt = Eigen::Vector3f::Zero();
    {
        Eigen::Vector3f t0 = vtx2xyz.row(i_vtx_fix) - vtx2xyz_ini.row(i_vtx_fix);
        dEdt += penalty * t0;
        dEdo += penalty * vtx2xyz_ini.row(i_vtx_fix).cross(t0.transpose() * rotation);
    }
    for (unsigned int i_vtx = 0; i_vtx < vtx2xyz.rows(); ++i_vtx)
    {
        // Write some code below to compute gradient of gravitational potential energy for each vertex
        // Code differentiation of energy w.r.t. translation and rotation for one line each.
        // For the differentiation w.r.t. rotation, observe how the rotation matrix will be updated at the line #83
        dEdt += -gravity;
        dEdo += -((rotation * vtx2xyz_ini.row(i_vtx).transpose()).cross(gravity));
        // do not change anything else except for the lines above.
    }
    translation -= learning_rate * dEdt;
    rotation = rotation * Eigen::AngleAxisf(-dEdo.norm() * learning_rate, dEdo.stableNormalized());
    vtx2xyz = ((rotation * vtx2xyz_ini.transpose()).colwise() + translation).transpose();
}

#define ACTIVATE_GAS_GEOMETRY static PRM_Name GeometryName(GAS_NAME_GEOMETRY, SIM_GEOMETRY_DATANAME); static PRM_Default GeometryNameDefault(0, SIM_GEOMETRY_DATANAME); PRMs.emplace_back(PRM_STRING, 1, &GeometryName, &GeometryNameDefault);
#define POINT_ATTRIBUTE_V3(NAME) GA_RWAttributeRef NAME##_attr = gdp.findGlobalAttribute(#NAME); if (!NAME##_attr.isValid()) NAME##_attr = gdp.addFloatTuple(GA_ATTRIB_POINT, #NAME, 3, GA_Defaults(0)); GA_RWHandleV3 NAME##_handle(NAME##_attr);

const SIM_DopDescription* GAS_Task09::getDopDescription()
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

void GAS_Task09::initializeSubclass()
{
    GAS_SubSolver::initializeSubclass();
    this->rotation = Eigen::Matrix3f::Identity();
    this->translation = Eigen::Vector3f::Zero();
}

void GAS_Task09::makeEqualSubclass(const SIM_Data* source)
{
    GAS_SubSolver::makeEqualSubclass(source);
    this->rotation = ((GAS_Task09 *) source)->rotation;
    this->translation = ((GAS_Task09 *) source)->translation;
}

bool GAS_Task09::solveGasSubclass(SIM_Engine& engine, SIM_Object* obj, SIM_Time time, SIM_Time timestep)
{
    SIM_GeometryCopy* G = getGeometryCopy(obj, GAS_NAME_GEOMETRY);
    SIM_GeometryAutoWriteLock lock(G);
    GU_Detail& gdp = lock.getGdp();
    POINT_ATTRIBUTE_V3(vtx2xyz_ini);

    Eigen::Matrix<float, Eigen::Dynamic, 3, Eigen::RowMajor> vtx2xyz;
    Eigen::Matrix<float, Eigen::Dynamic, 3, Eigen::RowMajor> vtx2xyz_ini;

    vtx2xyz.resize(gdp.getNumPoints(), 3);
    vtx2xyz_ini.resize(gdp.getNumPoints(), 3);

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

    constexpr unsigned int ivtx_fix = 0;
    step(rotation, translation, vtx2xyz, ivtx_fix, vtx2xyz_ini);

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
