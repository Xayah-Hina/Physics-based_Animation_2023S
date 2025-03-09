#ifndef TASK09_H
#define TASK09_H

#include <GAS/GAS_SubSolver.h>
#include <Eigen/Dense>


class GAS_Task09 final : public GAS_SubSolver
{
public:
    static constexpr bool GEN_NODE = true;
    inline static auto DOP_NAME = "Task09";
    inline static auto DOP_ENGLISH = "Task09";
    inline static auto DATANAME = "Task09";
    static constexpr bool UNIQUE_DATANAME = false;

    Eigen::Matrix3f rotation = Eigen::Matrix3f::Identity();
    Eigen::Vector3f translation = Eigen::Vector3f::Zero(); // translation to optimize

protected:
    explicit GAS_Task09(const SIM_DataFactory* factory): BaseClass(factory) {}
    bool solveGasSubclass(SIM_Engine& engine, SIM_Object* obj, SIM_Time time, SIM_Time timestep) override;
    static const SIM_DopDescription* getDopDescription();
    void initializeSubclass() override;
    void makeEqualSubclass(const SIM_Data* source) override;
    DECLARE_STANDARD_GETCASTTOTYPE();
    DECLARE_DATAFACTORY(GAS_Task09, GAS_SubSolver, "This is a Code for Task09", getDopDescription());
};



#endif //TASK09_H
