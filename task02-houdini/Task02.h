#ifndef TASK02_H
#define TASK02_H

#include <GAS/GAS_SubSolver.h>


class GAS_Task02 final : public GAS_SubSolver
{
public:
    static constexpr bool GEN_NODE = true;
    inline static auto DOP_NAME = "Task02";
    inline static auto DOP_ENGLISH = "Task02";
    inline static auto DATANAME = "Task02";
    static constexpr bool UNIQUE_DATANAME = false;

protected:
    explicit GAS_Task02(const SIM_DataFactory* factory): BaseClass(factory) {}
    bool solveGasSubclass(SIM_Engine& engine, SIM_Object* obj, SIM_Time time, SIM_Time timestep) override;
    static const SIM_DopDescription* getDopDescription();
    DECLARE_STANDARD_GETCASTTOTYPE();
    DECLARE_DATAFACTORY(GAS_Task02, GAS_SubSolver, "This is a Code for Task02", getDopDescription());
};



#endif //TASK02_H
