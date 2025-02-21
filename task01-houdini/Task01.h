#ifndef TASK01_H
#define TASK01_H

#include <GAS/GAS_SubSolver.h>


class GAS_Task01 final : public GAS_SubSolver
{
public:
    static constexpr bool GEN_NODE = true;
    inline static auto DOP_NAME = "Task01";
    inline static auto DOP_ENGLISH = "Task01";
    inline static auto DATANAME = "Task01";
    static constexpr bool UNIQUE_DATANAME = false;

protected:
    explicit GAS_Task01(const SIM_DataFactory* factory): BaseClass(factory) {}
    bool solveGasSubclass(SIM_Engine& engine, SIM_Object* obj, SIM_Time time, SIM_Time timestep) override;
    static const SIM_DopDescription* getDopDescription();
    DECLARE_STANDARD_GETCASTTOTYPE();
    DECLARE_DATAFACTORY(GAS_Task01, GAS_SubSolver, "This is a Code for Task01", getDopDescription());
};



#endif //TASK01_H
