#include "Task03.h"

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
    Eigen::Vector2f force; //! force from all the other particles
};

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

/**
 * position to grid coordinate
 * @param [in] pos input position
 * @param [in] box_size size of square box
 * @param [in] num_div number of division for grid
 * @return grid coordinate
 */
std::array<unsigned int, 2> pos2grid(
    const Eigen::Vector2f& pos,
    float box_size,
    unsigned int num_div)
{
    const auto fndiv = static_cast<float>(num_div);
    const auto ix = static_cast<unsigned int>((pos.x() / box_size + 0.5f) * fndiv);
    const auto iy = static_cast<unsigned int>((pos.y() / box_size + 0.5f) * fndiv);
    return {ix, iy};
}

/**
 * Gravitational force with softening
 * @param d relative position
 * @return force
 */
Eigen::Vector2f gravitational_force(const Eigen::Vector2f& d)
{
    constexpr float eps = 2.0e-3; // softening coefficient
    float r = sqrt(d.squaredNorm() + eps * eps);
    return (1.f / (r * r * r)) * d;
}

/**
 * For each particle, set summation of gravitational forces from all the other particles in a brute-force way O(N^2)
 * @param [in,out] particles
 */
void set_force_bruteforce(std::vector<Particle>& particles)
{
    for (unsigned int ip = 0; ip < particles.size(); ++ip)
    {
        particles[ip].force.setZero();
        for (unsigned int jp = 0; jp < particles.size(); ++jp)
        {
            if (ip == jp) { continue; }
            particles[ip].force += gravitational_force(particles[jp].pos - particles[ip].pos);
        }
    }
}

unsigned int abs_diff(unsigned int a, unsigned b)
{
    return a > b ? a - b : b - a;
}

// pair of particle index and grid index
struct ParticleGridIndex
{
    unsigned int particle_idx;
    unsigned int grid_idx;
};

// data structure for acceleration
class Acceleration
{
public:
    Acceleration(unsigned int num_particle, unsigned int num_div)
    {
        idx2pgi.resize(num_particle);
        grid2idx.resize(num_div * num_div + 1);
        grid2cg.resize(num_div * num_div);
    }

    std::vector<ParticleGridIndex> idx2pgi; // data of jagged array
    std::vector<unsigned int> grid2idx; // index of jagged array
    std::vector<Eigen::Vector2f> grid2cg; // the center of gravity of each grid
};

/**
 * For each particle, set summation of gravitational forces from all the other particles in an accelerated way
 * @param [in,out] particles particles
 * @param [in,out] acc acceleration data structure
 * @param [in] box_size box size
 * @param [in] num_div number of division of grid
 */
void set_force_accelerated(
    std::vector<Particle>& particles,
    Acceleration& acc,
    float box_size,
    unsigned int num_div)
{
    // computation of acceleration data structure
    acc.idx2pgi.clear();
    for (unsigned int ip = 0; ip < particles.size(); ++ip)
    {
        auto [ix, iy] = pos2grid(particles[ip].pos, box_size, num_div);
        acc.idx2pgi.push_back({ip, iy * num_div + ix});
    }
    std::sort(acc.idx2pgi.begin(), acc.idx2pgi.end(),
              [](const ParticleGridIndex& lhs, const ParticleGridIndex& rhs) { return lhs.grid_idx < rhs.grid_idx; });
    std::fill(acc.grid2idx.begin(), acc.grid2idx.end(), 0);
    // construct index of the jagged array
    for (auto& pg : acc.idx2pgi)
    {
        acc.grid2idx[pg.grid_idx + 1]++;
    }
    for (unsigned int i = 0; i < num_div * num_div; ++i)
    {
        acc.grid2idx[i + 1] += acc.grid2idx[i];
    }
    // compute the center of the gravity
    for (unsigned int ix = 0; ix < num_div; ++ix)
    {
        for (unsigned int iy = 0; iy < num_div; ++iy)
        {
            acc.grid2cg[iy * num_div + ix].setZero();
            for (unsigned int idx = acc.grid2idx[iy * num_div + ix]; idx < acc.grid2idx[iy * num_div + ix + 1]; ++idx)
            {
                unsigned int ip = acc.idx2pgi[idx].particle_idx;
                acc.grid2cg[iy * num_div + ix] += particles[ip].pos;
            }
            unsigned int np = acc.grid2idx[iy * num_div + ix + 1] - acc.grid2idx[iy * num_div + ix];
            if (np == 0) { continue; }
            acc.grid2cg[iy * num_div + ix] /= static_cast<float>(np);
        }
    }
    //
    for (unsigned int ip = 0; ip < particles.size(); ++ip)
    {
        auto [ix, iy] = pos2grid(particles[ip].pos, box_size, num_div); // grid coordinate of particle with index `ip`
        particles[ip].force.setZero(); // initialize force
        // loop over all the grid set force to the particle with index `ip`
        for (unsigned int jx = 0; jx < num_div; ++jx)
        {
            for (unsigned int jy = 0; jy < num_div; ++jy)
            {
                if (abs_diff(ix, jx) <= 1 && abs_diff(iy, jy) <= 1)
                {
                    // this grid is near to the particle as the grid indexes are close
                    for (unsigned int idx = acc.grid2idx[jy * num_div + jx]; idx < acc.grid2idx[jy * num_div + jx + 1]; ++idx)
                    {
                        assert(acc.idx2pgi[idx].grid_idx == jy * num_div + jx);
                        unsigned int jp = acc.idx2pgi[idx].particle_idx; // particle index
                        if (ip == jp) { continue; }
                        assert(pos2grid(particles[jp].pos, box_size, num_div)[0] == jx);
                        assert(pos2grid(particles[jp].pos, box_size, num_div)[1] == jy);
                        particles[ip].force += gravitational_force(particles[jp].pos - particles[ip].pos); // adding up forces from near particles
                    }
                }
                else
                {
                    // far field approximation
                    // write a few lines of code here to compute the force from far grid.
                    // use the center for the gravity of the grid : `acc.grid2cg[jy * num_div + jx]`
                }
            }
        }
    }
}


#define ACTIVATE_GAS_GEOMETRY static PRM_Name GeometryName(GAS_NAME_GEOMETRY, SIM_GEOMETRY_DATANAME); static PRM_Default GeometryNameDefault(0, SIM_GEOMETRY_DATANAME); PRMs.emplace_back(PRM_STRING, 1, &GeometryName, &GeometryNameDefault);
#define POINT_ATTRIBUTE_V3(NAME) GA_RWAttributeRef NAME##_attr = gdp.findGlobalAttribute(#NAME); if (!NAME##_attr.isValid()) NAME##_attr = gdp.addFloatTuple(GA_ATTRIB_POINT, #NAME, 3, GA_Defaults(0)); GA_RWHandleV3 NAME##_handle(NAME##_attr);
#define GLOBAL_ATTRIBUTE_F(NAME) GA_RWAttributeRef NAME##_attr = gdp.findGlobalAttribute(#NAME); if (!NAME##_attr.isValid()) NAME##_attr = gdp.addFloatTuple(GA_ATTRIB_DETAIL, #NAME, 1, GA_Defaults(0)); GA_RWHandleF NAME##_handle(NAME##_attr);

const SIM_DopDescription* GAS_Task03::getDopDescription()
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

bool GAS_Task03::solveGasSubclass(SIM_Engine& engine, SIM_Object* obj, SIM_Time time, SIM_Time timestep)
{
    SIM_GeometryCopy* G = getGeometryCopy(obj, GAS_NAME_GEOMETRY);
    SIM_GeometryAutoWriteLock lock(G);
    GU_Detail& gdp = lock.getGdp();
    POINT_ATTRIBUTE_V3(vel);
    POINT_ATTRIBUTE_V3(force);
    POINT_ATTRIBUTE_V3(Cd);

    const float box_size = 1.8;
    constexpr unsigned int num_div = 30;
    constexpr float dt = 0.00002f; // time step

    // particle information
    std::vector<Particle> particles(5000); // change the number of particles here
    if (engine.getSimulationFrame(time) == 1)
    {
        for (auto& p : particles)
        {
            // initialization
            p.pos.setRandom();
            p.pos *= box_size * 0.4f;
            // velocity is random but its magnitude is 1
            p.velo = {200 * p.pos.y(), -200 * p.pos.x()};
            // random RGB color
            p.color.setRandom();
            p.color = p.color * 0.5 + Eigen::Vector3f(0.5f, 0.5f, 0.5f);

            GA_Offset off = gdp.appendPoint();
            gdp.setPos3(off, UT_Vector3(p.pos.x(), 0, p.pos.y()));
            vel_handle.set(off, UT_Vector3(p.velo.x(), 0, p.velo.y()));
            force_handle.set(off, UT_Vector3(p.force.x(), 0, p.force.y()));
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
            particles[iter].force = Eigen::Vector2f(force_handle.get(off).x(), force_handle.get(off).z());
            particles[iter].color = Eigen::Vector3f(Cd_handle.get(off).x(), Cd_handle.get(off).y(), Cd_handle.get(off).z());
            ++iter;
        }
    }

    // switch brute-force/accelerated computation here by uncomment/comment below
    set_force_bruteforce(particles);
    // set_force_accelerated(particles, acceleration, box_size, num_div);

    for (auto& p : particles)
    {
        // leap frog time integration
        p.velo += p.force * dt; // update velocity
        p.pos += p.velo * dt; // update position
        // collision between particle and walls
        collision_circle_plane(p.pos, p.velo, 0.f, {-box_size * 0.5f, 0.f}, {+1.f, 0.f}); // left wall
        collision_circle_plane(p.pos, p.velo, 0.f, {0.f, -box_size * 0.5f}, {0.f, 1.f}); // bottom wall
        collision_circle_plane(p.pos, p.velo, 0.f, {+box_size * 0.5f, 0.f}, {-1.f, 0.f}); // right wall
        collision_circle_plane(p.pos, p.velo, 0.f, {0.f, +box_size * 0.5f}, {0.f, -1.f}); // top wall
    }

    {
        GA_Offset off;
        unsigned int iter = 0;
        GA_FOR_ALL_PTOFF(&gdp, off)
        {
            gdp.setPos3(off, UT_Vector3(particles[iter].pos.x(), 0, particles[iter].pos.y()));
            vel_handle.set(off, UT_Vector3(particles[iter].velo.x(), 0, particles[iter].velo.y()));
            force_handle.set(off, UT_Vector3(particles[iter].force.x(), 0, particles[iter].force.y()));
            Cd_handle.set(off, UT_Vector3(particles[iter].color.x(), particles[iter].color.y(), particles[iter].color.z()));
            ++iter;
        }
    }

    return true;
}
