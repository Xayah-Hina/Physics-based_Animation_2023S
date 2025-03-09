#include "Task07.h"

#include <SIM/SIM_Engine.h>
#include <SIM/SIM_Object.h>
#include <SIM/SIM_PositionSimple.h>
#include <SIM/SIM_GeometryCopy.h>
#include <SIM/SIM_DopDescription.h>
#include <SIM/SIM_ScalarField.h>
#include <PRM/PRM_Template.h>
#include <PRM/PRM_Default.h>
#include <PRM/PRM_Name.h>

#include <GU/GU_PrimPoly.h>
#include <Eigen/Dense>

#include <random>

void solve_laplace_gauss_seidel_on_grid(
    std::vector<float>& vtx2val,
    unsigned int grid_size,
    const std::vector<bool>& vtx2isfix)
{
    assert(vtx2val.size() == grid_size * grid_size);
    for (unsigned int ix = 0; ix < grid_size; ++ix)
    {
        for (unsigned int iy = 0; iy < grid_size; ++iy)
        {
            unsigned int idx_center = iy * grid_size + ix;
            if (vtx2isfix[idx_center]) { continue; }
            // write some code below to implement Gauss-Seidel method
            // Do not write more than 5 lines of code
            float v_left = vtx2val[idx_center - 1];
            float v_right = vtx2val[idx_center + 1];
            float v_up = vtx2val[idx_center + grid_size];
            float v_down = vtx2val[idx_center - grid_size];
            vtx2val[idx_center] = 0.25f * (v_left + v_right + v_up + v_down);
        }
    }
}

#define ACTIVATE_GAS_DENSITY static PRM_Name DensityName(GAS_NAME_DENSITY, "Density"); static PRM_Default DensityNameDefault(0, GAS_NAME_DENSITY); PRMs.emplace_back(PRM_STRING, 1, &DensityName, &DensityNameDefault);

const SIM_DopDescription* GAS_Task07::getDopDescription()
{
    static std::vector<PRM_Template> PRMs;
    PRMs.clear();
    ACTIVATE_GAS_DENSITY
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

template <typename T>
T TO_1D_IDX(const UT_Vector2T<T>& idx, const UT_Vector2T<T>& res) { return idx.x() + res.x() * idx.y(); }

void WriteField2DPartial(SIM_RawField* TARGET, const std::vector<float>& SOURCE, const int AXIS1, const int AXIS2, const UT_JobInfo& info)
{
    UT_VoxelArrayIteratorF vit;
    vit.setArray(TARGET->fieldNC());
    vit.setCompressOnExit(true);
    vit.setPartialRange(info.job(), info.numJobs());

    const UT_Vector2I res = {TARGET->getVoxelRes()[AXIS1], TARGET->getVoxelRes()[AXIS2]};
    for (vit.rewind(); !vit.atEnd(); vit.advance())
    {
        UT_Vector2I cell(vit.idx(AXIS1), vit.idx(AXIS2));
        const auto idx = TO_1D_IDX(cell, res);
        vit.setValue(static_cast<float>(SOURCE[idx]));
    }
}

THREADED_METHOD4(, true, WriteField2D, SIM_RawField *, TARGET, const std::vector<float>&, SOURCE, const int, AXIS1, const int, AXIS2);

bool GAS_Task07::solveGasSubclass(SIM_Engine& engine, SIM_Object* obj, SIM_Time time, SIM_Time timestep)
{
    SIM_ScalarField* D = getScalarField(obj, GAS_NAME_DENSITY);

    constexpr float box_size = 1.8;
    unsigned int grid_size = 256; // grid resolution
    std::vector<float> vtx2val(grid_size * grid_size); // grid data array storing distances
    std::vector<bool> vtx2isfix(grid_size * grid_size, false);

    if (engine.getSimulationFrame(time) == 1)
    {
        std::mt19937 rand(std::random_device{}());
        std::uniform_real_distribution<float> dist01(0.f, 1.f);
        for (float& val : vtx2val) { val = dist01(rand); }

        for (unsigned int i = 0; i < grid_size; ++i)
        {
            vtx2isfix[i] = true;
            vtx2isfix[i * grid_size] = true;
            vtx2isfix[i * grid_size + (grid_size - 1)] = true;
            vtx2isfix[(grid_size - 1) * grid_size + i] = true;
            vtx2val[i] = 0.0;
            vtx2val[i * grid_size] = 0.0;
            vtx2val[i * grid_size + (grid_size - 1)] = 0.0;
            vtx2val[(grid_size - 1) * grid_size + i] = 0.0;
        }
        float h = box_size / static_cast<float>(grid_size - 1);
        for (unsigned int ix = 0; ix < grid_size; ++ix)
        {
            for (unsigned int iy = 0; iy < grid_size; ++iy)
            {
                const float x = static_cast<float>(ix) * h - box_size * 0.5f;
                const float y = static_cast<float>(iy) * h - box_size * 0.5f;
                const float r = sqrt(x * x + y * y);
                if (r < box_size * 0.4 && r > box_size * 0.3)
                {
                    if (fabs(y) < box_size * 0.2 && x > 0.0)
                    {
                        continue;
                    }
                    vtx2isfix[iy * grid_size + ix] = true;
                    vtx2val[iy * grid_size + ix] = 0.9;
                }
            }
        }
    }
    else
    {
        for (unsigned int i = 0; i < grid_size; ++i)
        {
            vtx2isfix[i] = true;
            vtx2isfix[i * grid_size] = true;
            vtx2isfix[i * grid_size + (grid_size - 1)] = true;
            vtx2isfix[(grid_size - 1) * grid_size + i] = true;
        }
        float h = box_size / static_cast<float>(grid_size - 1);
        for (unsigned int ix = 0; ix < grid_size; ++ix)
        {
            for (unsigned int iy = 0; iy < grid_size; ++iy)
            {
                const float x = static_cast<float>(ix) * h - box_size * 0.5f;
                const float y = static_cast<float>(iy) * h - box_size * 0.5f;
                const float r = sqrt(x * x + y * y);
                if (r < box_size * 0.4 && r > box_size * 0.3)
                {
                    if (fabs(y) < box_size * 0.2 && x > 0.0)
                    {
                        continue;
                    }
                    vtx2isfix[iy * grid_size + ix] = true;
                }
            }
        }

        UT_VoxelArrayIteratorF vit;
        vit.setConstArray(D->getField()->field());
        const UT_Vector2I res = {D->getField()->getVoxelRes().x(), D->getField()->getVoxelRes().y()};
        for (vit.rewind(); !vit.atEnd(); vit.advance())
        {
            const UT_Vector2I cell(vit.x(), vit.y());
            int idx = static_cast<int>(TO_1D_IDX(cell, res));
            vtx2val[idx] = vit.getValue();
        }
    }
    solve_laplace_gauss_seidel_on_grid(vtx2val, grid_size, vtx2isfix);
    WriteField2D(D->getField(), vtx2val, 0, 1);

    return true;
}
