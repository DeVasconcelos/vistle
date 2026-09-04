#ifndef VISTLE_STREAMLINEVTKM_STREAMLINEVTKM_H
#define VISTLE_STREAMLINEVTKM_STREAMLINEVTKM_H

#include <mutex>
#include <vector>

#include <viskores/cont/DataSet.h>
#include <viskores/cont/PartitionedDataSet.h>
#include <viskores/filter/Filter.h>
#include <viskores/Particle.h>

#include <vistle/alg/objalg.h>
#include <vistle/module/module.h>
#include <vistle/util/enum.h>
#include <vistle/vtkm/module_status.h>

#include <viskores/cont/ArrayHandle.h>


// TODO: narrow conversion from vistle::Float to viskores::FloatDefault (or vistle::Float to Particles)
// TODO: is it possible to also do backwards integration (like in Tracer)?
//       --> from what I can tell, the filter only supports forward integration as neither the filter itself
//           nor the integrator worklets nor the Particle data structure itself have an option for backwards
//           integration. So we would have to implement that ourselves.

class StreamlineVtkm: public vistle::Module {
public:
    struct GlobalData {
        std::vector<viskores::cont::PartitionedDataSet> partitionedDatasets;
        std::vector<vistle::Object::const_ptr> outputGrids;
        std::vector<std::vector<vistle::DataBase::ptr>> outputFields;

        std::mutex mutex;

        void clear();
        bool isEmpty() const;
    };

    struct InputData {
        vistle::Object::const_ptr vistleGrid;
        std::vector<vistle::DataBase::const_ptr> fields;

        viskores::cont::DataSet viskoresDataset;
    };

    struct OutputData {
        vistle::Object::const_ptr vistleGrid;
        std::vector<vistle::DataBase::ptr> fields;

        viskores::cont::DataSet viskoresDataset;
    };

    DEFINE_ENUM_WITH_STRING_CONVERSIONS(MappedDataHandling, (Use)(Require)(Discard)(Generate))

    StreamlineVtkm(const std::string &name, int moduleID, mpi::communicator comm);

private:
    const int m_numPorts;
    const MappedDataHandling m_mappedDataHandling;

    std::vector<vistle::Port *> m_inputPorts, m_outputPorts;

    mutable GlobalData m_globalData;

    vistle::IntParameter *m_printObjectInfo = nullptr;

    std::string getFieldName(int index, bool output = false) const;

    bool prepare() override;

    bool reduce(int timestep) override;

    bool checkAndNotify(const ModuleStatusPtr &status) const;

    bool tryToExecuteFilter(viskores::filter::Filter &filter, const viskores::cont::DataSet &inputDataset,
                            viskores::cont::DataSet &outputDataset) const;

    bool tryToExecuteFilter(viskores::filter::Filter &filter, const viskores::cont::PartitionedDataSet &inputDataset,
                            viskores::cont::PartitionedDataSet &outputDataset) const;

    ModuleStatusPtr readInPorts(const std::shared_ptr<vistle::BlockTask> &task, vistle::Object::const_ptr &grid,
                                std::vector<vistle::DataBase::const_ptr> &fields) const;

    // ---------------------------------------------------------------------------------

    vistle::IntParameter *m_integrationMethod;
    vistle::IntParameter *m_numberOfSeeds, *m_maxNumberOfSeeds, *m_numberOfSteps;
    vistle::IntParameter *m_startStyle;

    vistle::FloatParameter *m_stepSize;

    vistle::VectorParameter *m_direction;
    vistle::VectorParameter *m_startPoint1, *m_startPoint2;

    bool changeParameter(const vistle::Parameter *param) override;

    ModuleStatusPtr prepareInputGrid(InputData &input) const;

    ModuleStatusPtr prepareInputField(const vistle::Port *port, InputData &input, int index) const;

    bool compute(const std::shared_ptr<vistle::BlockTask> &task) const override;

    std::unique_ptr<viskores::filter::Filter> setUpFilter() const;

    vistle::Object::const_ptr prepareOutputGrid(const InputData &input, OutputData &output) const;

    vistle::DataBase::ptr prepareOutputField(const InputData &input, OutputData &output, int index,
                                             const std::string &fieldName) const;

    viskores::cont::ArrayHandle<viskores::Particle> createSeedArray() const;
};

#endif // VISTLE_STREAMLINEVTKM_STREAMLINEVTKM_H
