#ifndef TASK06_H
#define TASK06_H

#include <GAS/GAS_SubSolver.h>


class GAS_Task06 final : public GAS_SubSolver
{
public:
    static constexpr bool GEN_NODE = true;
    inline static auto DOP_NAME = "Task06";
    inline static auto DOP_ENGLISH = "Task06";
    inline static auto DATANAME = "Task06";
    static constexpr bool UNIQUE_DATANAME = false;

    GETSET_DATA_FUNCS_I("SolverType", SolverType)

protected:
    explicit GAS_Task06(const SIM_DataFactory* factory): BaseClass(factory) {}
    bool solveGasSubclass(SIM_Engine& engine, SIM_Object* obj, SIM_Time time, SIM_Time timestep) override;
    static const SIM_DopDescription* getDopDescription();
    DECLARE_STANDARD_GETCASTTOTYPE();
    DECLARE_DATAFACTORY(GAS_Task06, GAS_SubSolver, "This is a Code for Task06", getDopDescription());
};



#endif //TASK06_H
