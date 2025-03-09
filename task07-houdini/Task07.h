#ifndef TASK07_H
#define TASK07_H

#include <GAS/GAS_SubSolver.h>


class GAS_Task07 final : public GAS_SubSolver
{
public:
    static constexpr bool GEN_NODE = true;
    inline static auto DOP_NAME = "Task07";
    inline static auto DOP_ENGLISH = "Task07";
    inline static auto DATANAME = "Task07";
    static constexpr bool UNIQUE_DATANAME = false;

protected:
    explicit GAS_Task07(const SIM_DataFactory* factory): BaseClass(factory) {}
    bool solveGasSubclass(SIM_Engine& engine, SIM_Object* obj, SIM_Time time, SIM_Time timestep) override;
    static const SIM_DopDescription* getDopDescription();
    DECLARE_STANDARD_GETCASTTOTYPE();
    DECLARE_DATAFACTORY(GAS_Task07, GAS_SubSolver, "This is a Code for Task07", getDopDescription());
};



#endif //TASK07_H
