#ifndef TASK08_H
#define TASK08_H

#include <GAS/GAS_SubSolver.h>


class GAS_Task08 final : public GAS_SubSolver
{
public:
    static constexpr bool GEN_NODE = true;
    inline static auto DOP_NAME = "Task08";
    inline static auto DOP_ENGLISH = "Task08";
    inline static auto DATANAME = "Task08";
    static constexpr bool UNIQUE_DATANAME = false;

protected:
    explicit GAS_Task08(const SIM_DataFactory* factory): BaseClass(factory) {}
    bool solveGasSubclass(SIM_Engine& engine, SIM_Object* obj, SIM_Time time, SIM_Time timestep) override;
    static const SIM_DopDescription* getDopDescription();
    DECLARE_STANDARD_GETCASTTOTYPE();
    DECLARE_DATAFACTORY(GAS_Task08, GAS_SubSolver, "This is a Code for Task08", getDopDescription());
};



#endif //TASK08_H
