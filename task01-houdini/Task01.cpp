#include "Task01.h"

#include <SIM/SIM_Object.h>
#include <SIM/SIM_GeometryCopy.h>
#include <SIM/SIM_DopDescription.h>
#include <PRM/PRM_Template.h>
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
    A << 1.f, 0.f, 0.f, 1.f;
    b << r0, v0;
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

const SIM_DopDescription* GAS_Task01::getDopDescription()
{
    static std::vector<PRM_Template> PRMs;
    PRMs.clear();
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

    UT_DMatrix4& xform = G->lockTransform();
    UT_Vector3 old;
    xform.getTranslates(old);
    xform.setTranslates(UT_Vector3(old.x(), old.y(), 0.0f));
    G->releaseTransform();

    return true;
}
