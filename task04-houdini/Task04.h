#ifndef TASK04_H
#define TASK04_H

#include <GAS/GAS_SubSolver.h>


class GAS_Task04 final : public GAS_SubSolver
{
public:
    static constexpr bool GEN_NODE = true;
    inline static auto DOP_NAME = "Task04";
    inline static auto DOP_ENGLISH = "Task04";
    inline static auto DATANAME = "Task04";
    static constexpr bool UNIQUE_DATANAME = false;

    GETSET_DATA_FUNCS_I("SolverType", SolverType)

protected:
    explicit GAS_Task04(const SIM_DataFactory* factory): BaseClass(factory) {}
    bool solveGasSubclass(SIM_Engine& engine, SIM_Object* obj, SIM_Time time, SIM_Time timestep) override;
    static const SIM_DopDescription* getDopDescription();
    DECLARE_STANDARD_GETCASTTOTYPE();
    DECLARE_DATAFACTORY(GAS_Task04, GAS_SubSolver, "This is a Code for Task04", getDopDescription());
};



#endif //TASK04_H
