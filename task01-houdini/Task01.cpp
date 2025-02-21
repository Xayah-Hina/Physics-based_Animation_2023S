#include "Task01.h"

#include <SIM/SIM_Engine.h>
#include <SIM/SIM_Object.h>
#include <SIM/SIM_GeometryCopy.h>
#include <SIM/SIM_DopDescription.h>
#include <PRM/PRM_Template.h>
#include <PRM/PRM_Default.h>
#include <PRM/PRM_Name.h>
#include <Eigen/Dense>

/**
 * explicit time integration
 * @param p0 input radius and its velocity
 * @param dt time step
 * @return output radius and its velocity
 */
Eigen::Vector2f time_integration_explicit(const Eigen::Vector2f& p0, float dt)
{
    float r0 = p0.x(); // current radius
    float v0 = p0.y(); // current radius velocity
    float f0 = -1.0f / (r0 * r0);
    return {r0 + dt * v0, v0 + dt * f0};
}

/**
 * implicit time integration (assignment)
 * @param p0 input radius and its velocity
 * @param dt time step
 * @return output radius and its velocity
 */
Eigen::Vector2f time_integration_implicit(const Eigen::Vector2f& p0, float dt)
{
    const float r0 = p0.x(); // current radius
    const float v0 = p0.y(); // current radius velocity
    const float dfdr = 2.f / (r0 * r0 * r0); // hint!
    float f0 = -1.0f / (r0 * r0); // force
    Eigen::Matrix2f A;
    Eigen::Vector2f b;
    // modify the following two lines to implement implicit time integration

    // 隐式欧拉积分方程:
    // r1 = r0 + dt * v1
    // v1 = v0 + dt * ( f0 + dfdr * (r1 - r0) )
    //
    // 整理后得到矩阵形式:
    // [ 1      -dt         ] [ r1 ]   = [ r0 ]
    // [ -dt*dfdr   1        ] [ v1 ]     [ v0 + dt*(f0 - dfdr*r0) ]

    A << 1.f,      -dt,
         -dt * dfdr, 1.f;
    b << r0,
         v0 + dt * (f0 - dfdr * r0);
    return A.inverse() * b;
}

/**
 * reflecting ball at the ground
 * @param p0 input radius and its velocity
 * @return output radius and its velocity
 */
Eigen::Vector2f reflection(const Eigen::Vector2f& p0)
{
    if (p0.x() > 0.5f) { return p0; }
    float r0 = p0.x();
    float v0 = p0.y();
    float E = 0.5f * v0 * v0 - 1.f / r0; // energy before reflection
    return {0.5, std::sqrt(2.f * E + 4.f)}; // energy preserving reflection
}

#define ACTIVATE_GAS_GEOMETRY static PRM_Name GeometryName(GAS_NAME_GEOMETRY, SIM_GEOMETRY_DATANAME); static PRM_Default GeometryNameDefault(0, SIM_GEOMETRY_DATANAME); PRMs.emplace_back(PRM_STRING, 1, &GeometryName, &GeometryNameDefault);
#define GLOBAL_ATTRIBUTE_F(NAME) GA_RWAttributeRef NAME##_attr = gdp.findGlobalAttribute(#NAME); if (!NAME##_attr.isValid()) NAME##_attr = gdp.addFloatTuple(GA_ATTRIB_DETAIL, #NAME, 1, GA_Defaults(0)); GA_RWHandleF NAME##_handle(NAME##_attr);

const SIM_DopDescription* GAS_Task01::getDopDescription()
{
    static std::vector<PRM_Template> PRMs;
    PRMs.clear();
    ACTIVATE_GAS_GEOMETRY
    static std::array<PRM_Name, 3> IntegratorType = {
        PRM_Name("0", "Explicit"),
        PRM_Name("1", "Implicit"),
        PRM_Name(nullptr),
    };
    static PRM_Name IntegratorTypeName("IntegratorType", "Integrator Type");
    static PRM_Default IntegratorTypeNameDefault(0);
    static PRM_ChoiceList CLIntegratorType(PRM_CHOICELIST_SINGLE, IntegratorType.data());
    PRMs.emplace_back(PRM_ORD, 1, &IntegratorTypeName, &IntegratorTypeNameDefault, &CLIntegratorType);
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

bool GAS_Task01::solveGasSubclass(SIM_Engine& engine, SIM_Object* obj, SIM_Time time, SIM_Time timestep)
{
    SIM_GeometryCopy* G = getGeometryCopy(obj, GAS_NAME_GEOMETRY);
    SIM_GeometryAutoWriteLock lock(G);
    GU_Detail& gdp = lock.getGdp();
    UT_Vector3 old_pos = gdp.getPos3(0);
    GLOBAL_ATTRIBUTE_F(vel);
    GLOBAL_ATTRIBUTE_F(time_shift);

    if (engine.getSimulationFrame(time) == 1)
        time_shift_handle.set(0, static_cast<float>(std::atan2(old_pos.z(), old_pos.x())));

    float dt = static_cast<float>(timestep);
    float time_shift = time_shift_handle.get(0);
    Eigen::Vector2f phase(
        static_cast<float>(old_pos.length()),
        static_cast<float>(vel_handle.get(0))
    );
    switch (getIntegratorType())
    {
    case 0:
        phase = time_integration_explicit(phase, dt);
        break;
    case 1:
        phase = time_integration_implicit(phase, dt);
        break;
    default:
        throw std::runtime_error("Invalid Integrator Type");
    }
    phase = reflection(phase);
    gdp.setPos3(0, UT_Vector3(phase.x() * std::cos(time + time_shift), old_pos.y(),
                              phase.x() * std::sin(time + time_shift)));
    vel_handle.set(0, phase.y());

    return true;
}
