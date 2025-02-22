#include "Task02.h"
#include "Task02.h"

#include <SIM/SIM_Engine.h>
#include <SIM/SIM_Object.h>
#include <SIM/SIM_PositionSimple.h>
#include <SIM/SIM_GeometryCopy.h>
#include <SIM/SIM_DopDescription.h>
#include <PRM/PRM_Template.h>
#include <PRM/PRM_Default.h>
#include <PRM/PRM_Name.h>
#include <Eigen/Dense>

/**
 * particle class (radius = 0)
 */
class Particle
{
public:
    Eigen::Vector2f pos; //! position
    Eigen::Vector2f velo; //! velocity
    Eigen::Vector3f color; //! color
};

/**
 * collision between a particle and a ball
 * @param [in,out] p particle class (position & velocity etc..)
 * @param [in] particle_mass particle mass
 * @param [in,out] ball_pos ball position
 * @param [in,out] ball_velo ball velocity
 * @param [in] ball_mass ball mass
 * @param [in] ball_rad ball rad
 */
void collide_particle_ball(
    Particle& p,
    float particle_mass,
    Eigen::Vector2f& ball_pos,
    Eigen::Vector2f& ball_velo,
    float ball_mass,
    float ball_rad)
{
    if ((p.pos - ball_pos).norm() > ball_rad) { return; }
    const Eigen::Vector2f plane_norm = (p.pos - ball_pos).normalized();
    const Eigen::Vector2f plane_org = ball_pos + plane_norm * ball_rad;
    float height = (p.pos - plane_org).dot(plane_norm);
    p.pos -= height * 2 * plane_norm;

    //--------------------------------------------------------------------- //
    // ASSIGNMENT: implement the collision between ball and particle.
    //             The coefficient of restitution is 1, e.g., the magnitude
    //             of relative velocity does not change after the collision.
    //             The linear momentum should be conserved.
    //             The force between ball and particle should only in the normal
    //             direction (e.g, `plane_norm`). In other words, there is
    //             no friction. You do not need to change the positions.

    // comment out the line below
    // p.velo -= 2.f * (p.velo - ball_velo).dot(plane_norm) * plane_norm;

    // write a few lines of code to compute the velocity of ball and particle
    // please uncomment the lines below
    const Eigen::Vector2f impulse = -(2.f * (p.velo - ball_velo).dot(plane_norm))
        / (1.f / particle_mass + 1.f / ball_mass) * plane_norm;
    p.velo += impulse / particle_mass;
    ball_velo -= impulse / ball_mass;
}

/**
 * collision between a circle and plane
 * @param [in,out] pos  circle position
 * @param [in, out] velo  circle velocity
 * @param [in] rad  circle radius
 * @param [in] plane_org  plane origin
 * @param [in] plane_nrm  plane normal
 */
void collision_circle_plane(
    Eigen::Vector2f& pos,
    Eigen::Vector2f& velo,
    const float rad,
    const Eigen::Vector2f& plane_org,
    const Eigen::Vector2f& plane_nrm)
{
    const float height = (pos - plane_org).dot(plane_nrm) - rad;
    if (height > 0.f) { return; }
    pos -= height * 2 * plane_nrm;
    const float velo_perp = velo.dot(plane_nrm);
    velo -= 2.f * velo_perp * plane_nrm;
}


#define ACTIVATE_GAS_GEOMETRY static PRM_Name GeometryName(GAS_NAME_GEOMETRY, SIM_GEOMETRY_DATANAME); static PRM_Default GeometryNameDefault(0, SIM_GEOMETRY_DATANAME); PRMs.emplace_back(PRM_STRING, 1, &GeometryName, &GeometryNameDefault);
#define ACTIVATE_GAS_GEOMETRY_BOX static PRM_Name GeometryBoxName("geo2", "GeometryBox"); static PRM_Default GeometryBoxNameDefault(0, "GeometryBox"); PRMs.emplace_back(PRM_STRING, 1, &GeometryBoxName, &GeometryBoxNameDefault);
#define ACTIVATE_GAS_GEOMETRY_CIRCLE static PRM_Name GeometryCircleName("geo3", "GeometryCircle"); static PRM_Default GeometryCircleNameDefault(0, "GeometryCircle"); PRMs.emplace_back(PRM_STRING, 1, &GeometryCircleName, &GeometryCircleNameDefault);
#define POINT_ATTRIBUTE_V3(NAME) GA_RWAttributeRef NAME##_attr = gdp.findGlobalAttribute(#NAME); if (!NAME##_attr.isValid()) NAME##_attr = gdp.addFloatTuple(GA_ATTRIB_POINT, #NAME, 3, GA_Defaults(0)); GA_RWHandleV3 NAME##_handle(NAME##_attr);
#define GLOBAL_ATTRIBUTE_F(NAME) GA_RWAttributeRef NAME##_attr = gdp.findGlobalAttribute(#NAME); if (!NAME##_attr.isValid()) NAME##_attr = gdp.addFloatTuple(GA_ATTRIB_DETAIL, #NAME, 1, GA_Defaults(0)); GA_RWHandleF NAME##_handle(NAME##_attr);

const SIM_DopDescription* GAS_Task02::getDopDescription()
{
    static std::vector<PRM_Template> PRMs;
    PRMs.clear();
    ACTIVATE_GAS_GEOMETRY
    ACTIVATE_GAS_GEOMETRY_BOX
    ACTIVATE_GAS_GEOMETRY_CIRCLE
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

bool GAS_Task02::solveGasSubclass(SIM_Engine& engine, SIM_Object* obj, SIM_Time time, SIM_Time timestep)
{
    SIM_GeometryCopy* G = getGeometryCopy(obj, GAS_NAME_GEOMETRY);
    SIM_GeometryCopy* BOX = getGeometryCopy(obj, "geo2");
    SIM_GeometryCopy* CIRCLE = getGeometryCopy(obj, "geo3");
    SIM_PositionSimple* POSITION_CIRCLE = SIM_DATA_GET(*obj, "PositionCircle", SIM_PositionSimple);


    const float box_size = 1.5;
    const float ball_rad = 0.3;
    const float ball_mass = 10.f;
    const float particle_mass = 1.0f;
    float dt = 0.01f;
    Eigen::Vector2f ball_pos(static_cast<float>(POSITION_CIRCLE->getPosition().x()), static_cast<float>(POSITION_CIRCLE->getPosition().z()));
    Eigen::Vector2f ball_velo;
    {
        SIM_GeometryAutoWriteLock lock_cir(CIRCLE);
        GU_Detail& gdp = lock_cir.getGdp();
        GLOBAL_ATTRIBUTE_F(Vx)
        GLOBAL_ATTRIBUTE_F(Vy)
        ball_velo = Eigen::Vector2f(static_cast<float>(Vx_handle.get(0)), static_cast<float>(Vy_handle.get(0)));
    }


    SIM_GeometryAutoWriteLock lock(G);
    GU_Detail& gdp = lock.getGdp();
    POINT_ATTRIBUTE_V3(vel);
    POINT_ATTRIBUTE_V3(Cd);

    std::vector<Particle> particles(100);
    if (engine.getSimulationFrame(time) == 1)
    {
        for (auto& p : particles)
        {
            // initialization
            // particle should be put outside the ball at the beginning
            for (unsigned int itr = 0; itr < 100; ++itr)
            {
                p.pos.setRandom();
                p.pos *= box_size * 0.5f;
                if (p.pos.norm() > ball_rad) { break; }
            }
            // velocity is random but its magnitude is 1
            p.velo.setRandom();
            p.velo.normalize();
            // random RGB color
            p.color.setRandom();
            p.color = p.color * 0.5 + Eigen::Vector3f(0.5f, 0.5f, 0.5f);

            GA_Offset off = gdp.appendPoint();
            gdp.setPos3(off, UT_Vector3(p.pos.x(), 0, p.pos.y()));
            vel_handle.set(off, UT_Vector3(p.velo.x(), 0, p.velo.y()));
            Cd_handle.set(off, UT_Vector3(p.color.x(), p.color.y(), p.color.z()));
        }
    }
    else
    {
        GA_Offset off;
        unsigned int iter = 0;
        GA_FOR_ALL_PTOFF(&gdp, off)
        {
            particles[iter].pos = Eigen::Vector2f(gdp.getPos3(off).x(), gdp.getPos3(off).z());
            particles[iter].velo = Eigen::Vector2f(vel_handle.get(off).x(), vel_handle.get(off).z());
            particles[iter].color = Eigen::Vector3f(Cd_handle.get(off).x(), Cd_handle.get(off).y(), Cd_handle.get(off).z());
            ++iter;
        }
    }

    ball_pos += ball_velo * dt; // step time for ball
    // collision between ball and the walls
    collision_circle_plane(ball_pos, ball_velo, ball_rad, {-box_size * 0.5f, 0.f}, {+1.f, 0.f}); // left wall
    collision_circle_plane(ball_pos, ball_velo, ball_rad, {0.f, -box_size * 0.5f}, {0.f, 1.f}); // bottom wall
    collision_circle_plane(ball_pos, ball_velo, ball_rad, {+box_size * 0.5f, 0.f}, {-1.f, 0.f}); // right wall
    collision_circle_plane(ball_pos, ball_velo, ball_rad, {0.f, +box_size * 0.5f}, {0.f, -1.f}); // top wall
    for (auto& p : particles)
    {
        p.pos += p.velo * dt; // step time for a particle
        // collision between a particle and the walls
        collision_circle_plane(p.pos, p.velo, 0.f, {-box_size * 0.5f, 0.f}, {+1.f, 0.f}); // left wall
        collision_circle_plane(p.pos, p.velo, 0.f, {0.f, -box_size * 0.5f}, {0.f, 1.f}); // bottom wall
        collision_circle_plane(p.pos, p.velo, 0.f, {+box_size * 0.5f, 0.f}, {-1.f, 0.f}); // right wall
        collision_circle_plane(p.pos, p.velo, 0.f, {0.f, +box_size * 0.5f}, {0.f, -1.f}); // top wall
        collide_particle_ball(p, particle_mass, ball_pos, ball_velo, ball_mass, ball_rad);
    }

    POSITION_CIRCLE->setPosition(UT_Vector3(ball_pos.x(), 0, ball_pos.y()));
    {
        SIM_GeometryAutoWriteLock lock_cir(CIRCLE);
        GU_Detail& gdp = lock_cir.getGdp();
        GLOBAL_ATTRIBUTE_F(Vx)
        GLOBAL_ATTRIBUTE_F(Vy)
        Vx_handle.set(0, ball_velo.x());
        Vy_handle.set(0, ball_velo.y());
    }
    {
        GA_Offset off;
        unsigned int iter = 0;
        GA_FOR_ALL_PTOFF(&gdp, off)
        {
            gdp.setPos3(off, UT_Vector3(particles[iter].pos.x(), 0, particles[iter].pos.y()));
            vel_handle.set(off, UT_Vector3(particles[iter].velo.x(), 0, particles[iter].velo.y()));
            Cd_handle.set(off, UT_Vector3(particles[iter].color.x(), particles[iter].color.y(), particles[iter].color.z()));
            ++iter;
        }
    }

    return true;
}
