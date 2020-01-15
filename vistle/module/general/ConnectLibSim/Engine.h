#ifndef VISIT_VISTLE_ENGINE_H
#define VISIT_VISTLE_ENGINE_H

#include <mpi.h>
#include <boost/asio.hpp>

#include "MetaData.h"
#include "VisItExports.h"

#include <string>
#include <vector>
#include <memory>
#include <functional>
#include <mutex>
#include <thread>

#ifdef MODULE_THREAD
#include <manager/manager.h>
#endif // MODULE_THREAD


#include <module/module.h>
#include <core/vec.h>

namespace in_situ {
enum class SimulationDataTyp {
     mesh
    ,variable
    ,material
    ,curve
    ,expression
    ,species
    ,genericCommand
    ,customCommand
    ,message

};

class V_VISITXPORT Engine {
public:
    static Engine* createEngine();
    static void DisconnectSimulation();
    bool initialize(int argC, char** argV);
    bool isInitialized() const noexcept;
    bool setMpiComm(void* newConn);


    //********************************
    //***functions called by module***
    //********************************
    void setModule(vistle::Module* module);
    void setDoReadMutex(std::mutex* m);
    int getNumObjects(SimulationDataTyp type);
    visit_handle getNthObject(SimulationDataTyp type, int n);
    std::vector<std::string> getDataNames(SimulationDataTyp type);
    //set callbacks (called from module)
    void SetTimestepChangedCb(std::function<bool(void)> cb);
    void SetDisconnectCb(std::function<void(void)> cb);
    //********************************
    //****functions called by sim****
    //********************************
    //adds all available data to the according outputs to execute the pipeline
    bool sendData();
    //called from simulation when a timestep changed
    void SimulationTimeStepChanged();
    void SimulationInitiateCommand(const char* command);
    void DeleteData();
    //set callbacks (called from sim)
    void SetSimulationCommandCallback(void(*sc)(const char*, const char*, void*), void* scdata);


private:
    static Engine* instance;
    bool m_initialized = false;
    MPI_Comm comm = MPI_COMM_WORLD;
    vistle::Module* m_module = nullptr;
    std::mutex* m_doReadMutex = nullptr;
    std::thread managerThread;

    std::map<std::string, vistle::Port*> m_portsList;
    struct MeshInfo {
        char* name = nullptr;
        int dim = 0; //2D or 3D
        int numDomains = 0;
        const int* domains = nullptr;
        std::vector<int> handles;
        std::vector< vistle::obj_ptr> grids;
    };
    std::map<std::string, MeshInfo> m_meshes;
    Metadata m_metaData;
    //callbacks from ConnectLibSim module
    std::function<bool(void)> timestepChangedCb; //returns true, if module is ready to receive data;
    std::function<void(void)> disconnectCb;


    //callbacks from simulation
    void (*simulationCommandCallback)(const char*, const char*, void*) = nullptr;
    void* simulationCommandCallbackData = nullptr;

    //retrieves metaData from simulation 
    void getMetaData();

    void addPorts();

    bool makeCurvilinearMesh(visit_handle meshMetaHandle);
    bool makeUntructuredMesh(visit_handle meshMetaHandle);
    bool makeAmrMesh(MeshInfo meshInfo);
    bool makeStructuredMesh(MeshInfo meshInfo);
    void sendMeshesToModule();
    void sendVarablesToModule();
    void sendDataToModule();
    void sendTestData();

    template<typename T>
    void sendVariableToModule(const std::string& name, vistle::obj_const_ptr mesh, int domain, const T* data, int size) {
        typename vistle::Vec<T, 1>::ptr variable(new typename vistle::Vec<T, 1>(size));
        memcpy(variable->x().data(), data, size * sizeof(T));
        variable->setGrid(mesh);
        //variable->setTimestep(timestep);
        variable->setBlock(domain);
        variable->setMapping(vistle::DataBase::Element);
        variable->addAttribute("_species", name);
        m_module->addObject(name, variable);
    }


    Engine();
    ~Engine();


    void printToConsole(const std::string& msg) const;
};





}

#endif // !VISIT_VISTLE_ENGINE_H
