#ifndef TASK10_H
#define TASK10_H

#include <GAS/GAS_SubSolver.h>
#include <Eigen/Dense>


class GAS_Task10 final : public GAS_SubSolver
{
public:
    static constexpr bool GEN_NODE = true;
    inline static auto DOP_NAME = "Task10";
    inline static auto DOP_ENGLISH = "Task10";
    inline static auto DATANAME = "Task10";
    static constexpr bool UNIQUE_DATANAME = false;

    Eigen::Vector3f Omega;
    Eigen::Matrix3f rotation = Eigen::Matrix3f::Identity();

protected:
    explicit GAS_Task10(const SIM_DataFactory* factory): BaseClass(factory) {}
    bool solveGasSubclass(SIM_Engine& engine, SIM_Object* obj, SIM_Time time, SIM_Time timestep) override;
    static const SIM_DopDescription* getDopDescription();
    void initializeSubclass() override;
    void makeEqualSubclass(const SIM_Data* source) override;
    DECLARE_STANDARD_GETCASTTOTYPE();
    DECLARE_DATAFACTORY(GAS_Task10, GAS_SubSolver, "This is a Code for Task10", getDopDescription());
};



#endif //TASK10_H
