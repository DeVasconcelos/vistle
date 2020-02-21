#include "Engine.h"

#include <mpi.h>

#include "VisItDataInterfaceRuntime.h"
#include "SimulationMetaData.h"
#include "VisItDataTypes.h"
#include "Exeption.h"
#include "TransformArray.h"


#include "MeshMetaData.h"
#include "VariableMetaData.h"
#include "MaterialMetaData.h"
#include "SpeciesMetaData.h"
#include "CurveMetaData.h"
#include "MessageMetaData.h"
#include "CommandMetaData.h"
#include "ExpressionMetaData.h"
#include "DomainList.h"
#include "RectilinearMesh.h"
#include "CurvilinearMesh.h"
#include "VariableData.h"

#include <boost/mpi.hpp>
#include <control/hub.h>


#include <module/module.h>

#include <core/rectilineargrid.h>
#include <core/structuredgrid.h>

#include <core/message.h>
#include <core/messagequeue.h>
#include <core/tcpmessage.h>

#include <util/sleep.h>
#include <util/directory.h>
#include <util/hostname.h>
#include <util/listenv4v6.h>


#include <ostream>


#ifdef  LIBSIM_DEBUG
#define DEBUG_CERR std::cerr << "Engine: " << " [" << m_rank << "/" << m_mpiSize << "] "
#else
struct DoNothing {
public:
    template<typename T>
    DoNothing& operator<<(const T& t) {
        return *this;
    }

    DoNothing& operator<<(std::ostream& (*pManip)(std::ostream&)) {
        return *this;
    }
};
static DoNothing doNothingInstance;
#define DEBUG_CERR doNothingInstance
#endif
#define CERR std::cerr << "Engine: " << " [" << m_rank << "/" << m_mpiSize << "] "


using std::string; using std::vector;
using std::endl;
using namespace insitu;
using namespace insitu::message;
namespace asio = boost::asio;
Engine* Engine::instance = nullptr;

Engine* Engine::createEngine() {
    if (!instance) {
        instance = new Engine{};
    }
    return instance;

}

int Engine::getNumObjects(SimulationDataTyp type) {

    std::function<int(visit_handle, int&)> getNum;
    switch (type) {
    case SimulationDataTyp::mesh:
    {
        getNum = simv2_SimulationMetaData_getNumMeshes;
    }
    break;
    case SimulationDataTyp::variable:
    {
        getNum = simv2_SimulationMetaData_getNumVariables;
    }
    break;
    case SimulationDataTyp::material:
    {
        getNum = simv2_SimulationMetaData_getNumMaterials;
    }
    break;
    case SimulationDataTyp::curve:
    {
        getNum = simv2_SimulationMetaData_getNumCurves;
    }
    break;
    case SimulationDataTyp::expression:
    {
        getNum = simv2_SimulationMetaData_getNumExpressions;
    }
    break;
    case SimulationDataTyp::species:
    {
        getNum = simv2_SimulationMetaData_getNumSpecies;
    }
    break;
    case SimulationDataTyp::genericCommand:
    {
        getNum = simv2_SimulationMetaData_getNumGenericCommands;
    }
    break;
    case SimulationDataTyp::customCommand:
    {
        getNum = simv2_SimulationMetaData_getNumCustomCommands;
    }
    break;
    case SimulationDataTyp::message:
    {
        getNum = simv2_SimulationMetaData_getNumMessages;
    }
    break;
    default:
        throw EngineExeption("getDataNames called with invalid type") << (int)type;

        break;
    }

    int num = 0;
    v2check(getNum, m_metaData.handle, num);
    return num;
}

visit_handle Engine::getNthObject(SimulationDataTyp type, int n) {
    std::function<int(visit_handle, int, visit_handle&)> getObj;
    switch (type) {
    case SimulationDataTyp::mesh:
    {
        getObj = simv2_SimulationMetaData_getMesh;
    }
    break;
    case SimulationDataTyp::variable:
    {
        getObj = simv2_SimulationMetaData_getVariable;
    }
    break;
    case SimulationDataTyp::material:
    {
        getObj = simv2_SimulationMetaData_getMaterial;
    }
    break;
    case SimulationDataTyp::curve:
    {
        getObj = simv2_SimulationMetaData_getCurve;
    }
    break;
    case SimulationDataTyp::expression:
    {
        getObj = simv2_SimulationMetaData_getExpression;
    }
    break;
    case SimulationDataTyp::species:
    {
        getObj = simv2_SimulationMetaData_getSpecies;
    }
    break;
    case SimulationDataTyp::genericCommand:
    {
        getObj = simv2_SimulationMetaData_getGenericCommand;
    }
    break;
    case SimulationDataTyp::customCommand:
    {
        getObj = simv2_SimulationMetaData_getCustomCommand;
    }
    break;
    case SimulationDataTyp::message:
    {
        getObj = simv2_SimulationMetaData_getMessage;
    }
    break;
    default:
        throw EngineExeption("getDataNames called with invalid type");
        break;
    }
    visit_handle obj = m_metaData.handle;
    v2check(getObj, obj, n, obj);
    return obj;
}

std::vector<std::string> insitu::Engine::getDataNames(SimulationDataTyp type) {
    std::function<int(visit_handle, char**)> getName;
    switch (type) {
    case SimulationDataTyp::mesh:
    {
        getName = simv2_MeshMetaData_getName;
    }
    break;
    case SimulationDataTyp::variable:
    {
        getName = simv2_VariableMetaData_getName;
    }
    break;
    case SimulationDataTyp::material:
    {
        getName = simv2_MaterialMetaData_getName;
    }
    break;
    case SimulationDataTyp::curve:
    {
        getName = simv2_CurveMetaData_getName;
    }
    break;
    case SimulationDataTyp::expression:
    {
        getName = simv2_ExpressionMetaData_getName;
    }
    break;
    case SimulationDataTyp::species:
    {
        getName = simv2_SpeciesMetaData_getName;
    }
    break;
    case SimulationDataTyp::genericCommand:
    {
        getName = simv2_CommandMetaData_getName;
    }
    break;
    case SimulationDataTyp::customCommand:
    {
        getName = simv2_CommandMetaData_getName;
    }
    break;
    case SimulationDataTyp::message:
    {
        getName = simv2_MessageMetaData_getName;
    }
    break;
    default:
        throw new EngineExeption("getDataNames called with invalid type");
        break;
    }

    int n = getNumObjects(type);
    std::vector<string> names;
    names.reserve(n);
    for (size_t i = 0; i < n; i++) {
        visit_handle obj = getNthObject(type, i);
        char* name;
        v2check(getName, obj, &name);
        names.push_back(name);
    }
    return names;
}

void Engine::DisconnectSimulation() {

    delete instance;
    instance = nullptr;
}

bool Engine::initialize(int argC, char** argV) {
    MPI_Comm_rank(MPI_COMM_WORLD, &m_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &m_mpiSize);
    CERR<< "__________Engine args__________" << endl;
    for (size_t i = 0; i < argC; i++) {
        CERR << argV[i] << endl;
    }
    CERR << "_______________________________" << endl;
#ifdef MODULE_THREAD
    // start manager on cluster
    const char *VISTLE_ROOT = getenv("VISTLE_ROOT");
    if (!VISTLE_ROOT) {
        CERR << "VISTLE_ROOT not set to the path of the Vistle build directory." << endl;
        return false;
    }

    m_managerThread = std::thread([argC, argV, VISTLE_ROOT]() {
        std::string cmd{VISTLE_ROOT};
        cmd += "/bin/vistle_manager";
        std::vector<char*> args;
        args.push_back(const_cast<char*>(cmd.c_str()));
        for (int i = 1; i < argC; ++i) {
            args.push_back(argV[i]);
        }
        vistle::VistleManager manager;
        manager.run(args.size(), args.data());
    });
    m_initialized = true;

#else

    if (argC != 6) {
        CERR << "simulation requires exactly 6 parameters" << endl;
        return false;
    }

    if (atoi(argV[0]) != m_mpiSize) {
        CERR << "mpi size of simulation must match vistle's mpi size" << endl;
        return false;
    }
    m_shmName = argV[1];
    m_moduleName = argV[2];
    m_moduleID = atoi(argV[3]);
    if (m_rank == 0 && argV[4] != vistle::hostname())         {
        CERR << "this " << vistle::hostname() << "trying to connect to " << argV[4] << endl;
        CERR << "Wrong host: must connect to Vistle on the same machine!" << endl;
        return false;
    }

    vistle::registerTypes();
    vistle::Shm::attach(m_shmName, m_moduleID, m_rank);

#ifndef MODULE_THREAD
    vistle::message::DefaultSender::init(m_moduleID, m_rank);
#endif
    // names are swapped relative to communicator
    std::string mqName = vistle::message::MessageQueue::createName("recvFromSim", m_moduleID, m_rank);
    try {
        m_sendMessageQueue = vistle::message::MessageQueue::open(mqName);
    } catch (boost::interprocess::interprocess_exception & ex) {
        CERR << "opening send message queue " << mqName << ": " << ex.what() << endl;
       return false;
    }



    m_initialized = true;
#endif
    if (m_rank == 0) {
        try {
            initializeEngineSocket(argV[4], atoi(argV[5]));
        } catch (const EngineExeption& ex) {
            CERR << ex.what() << endl;
            return false;
        }
    }
    InSituTcpMessage::initialize(m_socket, boost::mpi::communicator(comm, boost::mpi::comm_create_kind::comm_duplicate));
    InSituTcpMessage::send(GoOn());
    SyncShmMessage::initialize(m_moduleID, m_rank, SyncShmMessage::Mode::Attach);

    return true;
}

bool Engine::isInitialized() const noexcept {
    return m_initialized;
}

bool Engine::setMpiComm(void* newconn) {
    comm = (MPI_Comm)newconn;

return true;

}

bool Engine::sendData() {

    CERR << "sendData was called" << endl;
    return true;
}

void Engine::SimulationTimeStepChanged() {
    if (!m_initialized || !m_moduleInitialized) {
        CERR << "not connected with Vistle. \nStart the ConnectLibSIm module to use simulation data in Vistle!" << endl;
        return;
    }
    static int counter = 0;
    DEBUG_CERR << "SimulationTimeStepChanged counter = " << counter << endl;
    ++counter;
    int oldCycle = m_metaData.currentCycle;
    getMetaData();
    if (oldCycle == m_metaData.currentCycle) {
        CERR << "There is no new timestep but SimulationTimeStepChanged was called" << endl;
        return;
    }
    int numMeshes, numVars;
    v2check(simv2_SimulationMetaData_getNumMeshes, m_metaData.handle, numMeshes);
    v2check(simv2_SimulationMetaData_getNumVariables, m_metaData.handle, numVars);
    if (m_metaData.currentCycle % m_nthTimestep != 0) {
        return;
    }
    DEBUG_CERR << "Timestep " << m_metaData.currentCycle << " has " << numMeshes << " meshes and " << numVars << "variables" << endl;
    if (m_moduleReady) {//only here vistle::objects are allowed to be made
        sendDataToModule();
    }
    else {
        CERR << "ConnectLibSim is not ready to process data" << endl;
    }
}

void Engine::SimulationInitiateCommand(const std::string& command) {
    static int counter = 0;
    DEBUG_CERR << "SimulationInitiateCommand " << command << " counter = " << counter << endl;
    ++counter;
    if (command.substr(0, 12) == "INTERNALSYNC") { //need to respond or LibSim gets stuck. see: SimEngine::SimulationInitiateCommand
        // Send the command back to the engine so it knows we're done syncing. 
        std::string cmd("INTERNALSYNC");
        std::string args(command.substr(13, command.size() - 1));
        simulationCommandCallback(cmd.c_str(), args.c_str(), simulationCommandCallbackData);
        InSituTcpMessage::send(GoOn{}); //request tcp message from conroller
    }
}

void Engine::DeleteData() {
    m_meshes.clear();
    sendData();
}

bool insitu::Engine::handleVistleMessage() {
    getMetaData();
    if (m_rank == 0) {
        if (!slaveCommandCallback) {
            CERR << "passCommandToSim failed : slaveCommandCallback not set" << endl;
            return true;;
        }
        slaveCommandCallback(); //let the slaves call visitProcessEngineCommand() -> simv2_process_input() -> Engine::passCommandToSim() and therefore finalizeInit() if not already done
    }
    static int counter = 0;
    DEBUG_CERR << "handleVistleMessage counter = " << counter << endl;
    ++counter;
    finalizeInit();
    InSituTcpMessage msg = InSituTcpMessage::recv();
    DEBUG_CERR << "received message of type " << static_cast<int>(msg.type()) << endl;
    switch (msg.type()) {
    case InSituMessageType::Invalid:
        break;
    case InSituMessageType::ShmInit:
        break;
    case InSituMessageType::AddObject:
        break;
    case InSituMessageType::AddPorts:
        break;
    case InSituMessageType::AddCommands:
        break;
    case InSituMessageType::Ready:
    {
        Ready em = msg.unpackOrCast<Ready>();
        m_moduleReady = em.m_state;
        if (m_moduleReady) {
            SyncShmMessage msg = SyncShmMessage::recv();
            vistle::Shm::the().setObjectID(msg.objectID());
            vistle::Shm::the().setArrayID(msg.arrayID());
        }
        else {
            SyncShmMessage::send(SyncShmMessage{ vistle::Shm::the().objectID(), vistle::Shm::the().arrayID() });
        }
    }
    break;
    case InSituMessageType::ExecuteCommand:
    {
        ExecuteCommand exe = msg.unpackOrCast<ExecuteCommand>();
        if (simulationCommandCallback) {
            simulationCommandCallback(exe.m_command.c_str(), "", simulationCommandCallbackData);
            DEBUG_CERR << "received simulation command: " << exe.m_command << endl;
            if (m_registeredGenericCommands.find(exe.m_command) == m_registeredGenericCommands.end()) {
                DEBUG_CERR << "Engine received unknown command!" << endl;
            }
        } else {

            CERR << "received command, but required callback is not set" << endl;
        }
        InSituTcpMessage::send(GoOn{});
    }
    break;
    case InSituMessageType::GoOn:
    {
        if (m_metaData.simMode == VISIT_SIMMODE_RUNNING) {
            //InSituTcpMessage::send(GoOn{});
        }
    }
        break;
    case InSituMessageType::ConstGrids:
    {
        ConstGrids em = msg.unpackOrCast<ConstGrids>();
        m_constGrids = em.m_state;
    }
    break;
    case InSituMessageType::NthTimestep:
    {
        NthTimestep em = msg.unpackOrCast<NthTimestep>();
        m_constGrids = em.m_frequency;
    }
    break;
    case InSituMessageType::ConnectionClosed:
    {
        CERR << "connection closed" << endl;
        return false;
    }
    break;
    default:
        break;
    }
    return true;
}

void Engine::SetSimulationCommandCallback(void(*sc)(const char*, const char*, void*), void* scdata) {
    simulationCommandCallback = sc;
    simulationCommandCallbackData = scdata;
}

void insitu::Engine::setSlaveComandCallback(void(*sc)(void)) {
    DEBUG_CERR << "setSlaveComandCallback" << endl;
    slaveCommandCallback = sc;
}

int insitu::Engine::GetInputSocket() {
    if (m_rank == 0 && m_socket) {
        return m_socket->native_handle();
    }
    else {
        return 0;
    }
}

void insitu::Engine::getMetaData() {
    m_metaData.handle = simv2_invoke_GetMetaData(); //somehow 0 is valid
    v2check(simv2_SimulationMetaData_getData, m_metaData.handle, m_metaData.simMode, m_metaData.currentCycle, m_metaData.currentTime);
}

void insitu::Engine::getRegisteredGenericCommands() {
    int numRegisteredCommands = getNumObjects(SimulationDataTyp::genericCommand);
    bool found = false;
    std::vector<std::string> commands;
    for (size_t i = 0; i < numRegisteredCommands; i++) {
        visit_handle commandHandle = getNthObject(SimulationDataTyp::genericCommand, i);
        char* name;
        v2check(simv2_CommandMetaData_getName, commandHandle, &name);
        CERR << "registerd generic command: " << name << endl;
        m_registeredGenericCommands.insert(name);
        commands.push_back(name);
    }
    InSituTcpMessage::send(AddCommands{ commands });
}

void Engine::addPorts() {
    try {
        std::vector<string> names = getDataNames(SimulationDataTyp::mesh);
        names.push_back("mesh");
        InSituTcpMessage::send(AddPorts{ names });

        names = getDataNames(SimulationDataTyp::variable);
        names.push_back("variable");
        InSituTcpMessage::send(AddPorts{ names });

    } catch (const SimV2Exeption& ex) {
        CERR << "failed to add output ports: " << ex.what() << endl;
    } catch (const EngineExeption & ex) {
        CERR << "failed to add output ports: " << ex.what() << endl;
    }
}

bool insitu::Engine::makeRectilinearMesh(MeshInfo meshInfo) {
    for (size_t cd = 0; cd < meshInfo.numDomains; cd++) {
        int currDomain = meshInfo.domains[cd];
        visit_handle meshHandle = v2check(simv2_invoke_GetMesh, currDomain, meshInfo.name);
        int check = simv2_RectilinearMesh_check(meshHandle);
        if (check == VISIT_OKAY) {
            visit_handle coordHandles[3]; //handles to variable data
            int ndims;
            v2check(simv2_RectilinearMesh_getCoords, meshHandle, &ndims, &coordHandles[0], &coordHandles[1], &coordHandles[2]);
            std::array<int, 3> owner{}, dataType{}, nComps{}, nTuples{ 1,1,1 };
            std::array<void*, 3> data{};
            for (int i = 0; i < meshInfo.dim; ++i) {
                v2check(simv2_VariableData_getData, coordHandles[i], owner[i], dataType[i], nComps[i], nTuples[i], data[i]);
                if (dataType[i] != dataType[0]) {
                    CERR << "mesh data type must be consisten within a domain" << endl;
                    return false;
                }
            }
            std::reverse(owner.begin(), owner.end());
            std::reverse(dataType.begin(), dataType.end());
            std::reverse(nComps.begin(), nComps.end());
            std::reverse(nTuples.begin(), nTuples.end());
            std::reverse(data.begin(), data.end());

            vistle::RectilinearGrid::ptr grid = vistle::RectilinearGrid::ptr(new vistle::RectilinearGrid(nTuples[0], nTuples[1], nTuples[2]));
            grid->setTimestep(m_constGrids ? -1 : m_metaData.currentCycle / m_nthTimestep);
            grid->setBlock(currDomain);
            meshInfo.handles.push_back(meshHandle);
            meshInfo.grids.push_back(grid);
            for (size_t i = 0; i < 3; ++i) {
                if (data[i]) {
                    transformArray(data[i], grid->coords(i).begin(), nTuples[i], dataType[i]);
                }
                else {
                    grid->coords(i)[0] = 0;
                }
            }

            std::array<int, 3> min, max;
            v2check(simv2_RectilinearMesh_getRealIndices, meshHandle, min.data(), max.data());
            std::reverse(min.begin(), min.end());
            std::reverse(max.begin(), max.end());
            for (size_t i = 0; i < 3; i++) {
                assert(min[i] >= 0);
                int numTop = nTuples[i] - 1 - max[i];
                assert(numTop >= 0);
                grid->setNumGhostLayers(i, vistle::StructuredGrid::GhostLayerPosition::Bottom, min[i]);
                grid->setNumGhostLayers(i, vistle::StructuredGrid::GhostLayerPosition::Top, numTop);
            }
            addObject(meshInfo.name, grid);

            DEBUG_CERR << "added rectilinear mesh " << meshInfo.name << " dom = " << currDomain << " xDim = " << nTuples[0] << " yDim = " << nTuples[1] << " zDim = " << nTuples[2] << endl;
            DEBUG_CERR << "min bounds " << min[0] << " " << min[1] << " " << min[2] << " max bouns " << max[0] << " " << max[1] << " " << max[2] << endl;

        }
    }
    DEBUG_CERR << "made rectilinear grids with " << meshInfo.numDomains << " domains" << endl;
    m_meshes[meshInfo.name] = meshInfo;
    return true;


}

bool insitu::Engine::makeUntructuredMesh(MeshInfo meshInfo) {

    return false;
}

bool insitu::Engine::makeAmrMesh(MeshInfo meshInfo) {
    return makeRectilinearMesh(meshInfo);

}

bool insitu::Engine::makeStructuredMesh(MeshInfo meshInfo) {
    for (size_t cd = 0; cd < meshInfo.numDomains; cd++) {
        int currDomain = meshInfo.domains[cd];
        visit_handle meshHandle = v2check(simv2_invoke_GetMesh, currDomain, meshInfo.name);
        int check = simv2_CurvilinearMesh_check(meshHandle);
        if (check == VISIT_OKAY) {
            visit_handle coordHandles[4]; //handles to variable data, 4th entry conteins interleaved data depending on coordMode
            int dims[3]{ 1,1,1 }; //the x,y,z dimensions
            int ndims, coordMode;
            //no v2check because last visit_handle can be invalid
            if (!simv2_CurvilinearMesh_getCoords(meshHandle, &ndims, dims, &coordMode, &coordHandles[0], &coordHandles[1], &coordHandles[2], &coordHandles[3]))                 {
                throw EngineExeption("makeStructuredMesh: simv2_CurvilinearMesh_getCoords failed");
            }
            vistle::StructuredGrid::ptr grid(new vistle::StructuredGrid(dims[0], dims[1], dims[2]));
            std::array<float*, 3> gridCoords{ grid->x().data() ,grid->y().data() ,grid->z().data() };
            int numVals = dims[0] * dims[1] * dims[2];
            int owner{}, dataType{}, nComps{}, nTuples{};
            void* data{};
            switch (coordMode) {
            case VISIT_COORD_MODE_INTERLEAVED:
            {
                v2check(simv2_VariableData_getData, coordHandles[3], owner, dataType, nComps, nTuples, data);
                if (nTuples != numVals * ndims) {
                    throw EngineExeption("makeStructuredMesh: received points in interleaved grid " + std::string(meshInfo.name) + " do not match the grid dimensions");
                }
               
                transformInterleavedArray(data, gridCoords, numVals, dataType, ndims);

            }
            break;
            case VISIT_COORD_MODE_SEPARATE:
            {
            for (int i = 0; i < meshInfo.dim; ++i) {
                v2check(simv2_VariableData_getData, coordHandles[i], owner, dataType, nComps, nTuples, data);
                if (nTuples != numVals) {
                    throw EngineExeption("makeStructuredMesh: received points in separate grid " + std::string(meshInfo.name) + " do not match the grid dimension " + std::to_string(i));
                }
                transformArray(data, gridCoords[i], nTuples, dataType);
            }
            }
            break;
            default:
                throw EngineExeption("coord mode must be interleaved(1) or separate(0), it is " + std::to_string(coordMode));
            }
            if (meshInfo.dim == 2) {
                std::fill(gridCoords[2], gridCoords[2] + numVals, 0);
            }
            grid->setTimestep(m_constGrids ? -1 : m_metaData.currentCycle / m_nthTimestep);
            grid->setBlock(currDomain);
            meshInfo.handles.push_back(meshHandle);
            meshInfo.grids.push_back(grid);
            addObject(meshInfo.name, grid);
        }
    }
    CERR << "made structured mesh with " << meshInfo.numDomains << " domains" << endl;
    m_meshes[meshInfo.name] = meshInfo;
    return true;
}

void insitu::Engine::sendMeshesToModule()     {
    int numMeshes = getNumObjects(SimulationDataTyp::mesh);

    for (size_t i = 0; i < numMeshes; i++) {
        MeshInfo meshInfo;
        visit_handle meshHandle = getNthObject(SimulationDataTyp::mesh, i);
        char* name;
        v2check(simv2_MeshMetaData_getName, meshHandle, &name);
        if (m_constGrids && m_meshes.find(name) != m_meshes.end()) {
            auto m = m_meshes.find(name);
            if (m != m_meshes.end()) {
                for (auto grid : m->second.grids)                     {
                    addObject(name, grid);
                }
                return; //the mesh is consistent and already there
            }
        }
        visit_handle domainListHandle = v2check(simv2_invoke_GetDomainList, name);
        int allDoms = 0;
        visit_handle myDoms;
        v2check(simv2_DomainList_getData, domainListHandle, allDoms, myDoms);
        DEBUG_CERR << "allDoms = " << allDoms << " myDoms = " << myDoms << endl;

        int owner, dataType, nComps;
        void* data;
        try {
            v2check(simv2_VariableData_getData, myDoms, owner, dataType, nComps, meshInfo.numDomains, data);
        } catch (const SimV2Exeption&) {
            CERR << "failed to get domain list" << endl;
        }
        if (dataType != VISIT_DATATYPE_INT) {
            CERR << "expected domain list to be ints" << endl;
        }
        meshInfo.domains = static_cast<const int*>(data);
  

        VisIt_MeshType meshType = VisIt_MeshType::VISIT_MESHTYPE_UNKNOWN;
        v2check(simv2_MeshMetaData_getMeshType, meshHandle, (int*)&meshType);


        v2check(simv2_MeshMetaData_getName, meshHandle, &meshInfo.name);
        v2check(simv2_MeshMetaData_getTopologicalDimension, meshHandle, &meshInfo.dim);

        switch (meshType) {
        case VISIT_MESHTYPE_RECTILINEAR:
        {
            CERR << "making rectilinear grid" << endl;
            makeRectilinearMesh(meshInfo);
        }
        break;
        case VISIT_MESHTYPE_CURVILINEAR:
        {
            CERR << "making curvilinear grid" << endl;
            makeStructuredMesh(meshInfo);
        }
        break;
        case VISIT_MESHTYPE_UNSTRUCTURED:
        {
            CERR << "making unstructured grid" << endl;
        }
        break;
        case VISIT_MESHTYPE_POINT:
        {
            CERR << "making point grid" << endl;
        }
        break;
        case VISIT_MESHTYPE_CSG:
        {
            CERR << "making csg grid" << endl;
        }
        break;
        case VISIT_MESHTYPE_AMR:
        {
            makeAmrMesh(meshInfo);
        }
        break;
        default:
            throw EngineExeption("unknown meshtype");
            break;
        }
    }
}

void insitu::Engine::sendVarablesToModule()     { //todo: combine variables to vectors
    int numVars = getNumObjects(SimulationDataTyp::variable);
    for (size_t i = 0; i < numVars; i++) {
        visit_handle varMetaHandle = getNthObject(SimulationDataTyp::variable, i);
        char* name, *meshName;
        v2check(simv2_VariableMetaData_getName, varMetaHandle, &name);
        v2check(simv2_VariableMetaData_getMeshName, varMetaHandle, &meshName);
        int centering = -1;
        v2check(simv2_VariableMetaData_getCentering, varMetaHandle, &centering);
        auto meshInfo = m_meshes.find(meshName);
        if (meshInfo == m_meshes.end()) {
            throw EngineExeption(std::string("can't find mesh ") + meshName + " for variable " + name);
        }
        for (size_t cd = 0; cd < meshInfo->second.numDomains; ++cd) {
            int currDomain = meshInfo->second.domains[cd];
            visit_handle varHandle = v2check(simv2_invoke_GetVariable, currDomain, name);
            int  owner{}, dataType{}, nComps{}, nTuples{};
            void* data = nullptr;
            v2check(simv2_VariableData_getData, varHandle, owner, dataType, nComps, nTuples, data);
            //CERR << "variable " << name << " domain " << currDomain << " owner = " << owner << " dataType = " << dataType << " ncomps = " << nComps << " nTuples = " << nTuples  << endl;
            vistle::Vec<vistle::Scalar, 1>::ptr variable(new typename vistle::Vec<vistle::Scalar, 1>(nTuples));
            transformArray(data, variable->x().data(), nTuples, dataType);
            variable->setGrid(meshInfo->second.grids[cd]);
            variable->setTimestep(m_metaData.currentCycle / m_nthTimestep);
            variable->setBlock(currDomain);
            variable->setMapping(centering == VISIT_VARCENTERING_NODE? vistle::DataBase::Vertex : vistle::DataBase::Element);
            variable->addAttribute("_species", name);
            addObject(name, variable);
            DEBUG_CERR << "added variable " << name  << " to mesh " << meshName << " dom = " << currDomain << " Dim = " << nTuples << endl;

        }
        DEBUG_CERR << "sent variable " << name << " with " << meshInfo->second.numDomains << " domains" << endl;
    }
}

void insitu::Engine::sendDataToModule() {
    try {
        sendMeshesToModule();
    } catch (const VistleLibSimExeption & exept) {
        CERR << "sendMeshesToModule failed: " << exept.what() << endl;
    }

    try {
        sendVarablesToModule();
    } catch (const VistleLibSimExeption & exept) {
        CERR << "sendVarablesToModule failed: " << exept.what() << endl;
    }
}

void insitu::Engine::sendTestData() {

    int rank = -1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::array<int, 3> dims{ 5, 10, 1 };
    vistle::RectilinearGrid::ptr grid = vistle::RectilinearGrid::ptr(new vistle::RectilinearGrid(dims[0], dims[1], dims[2]));
    grid->setTimestep(m_constGrids ? -1 : m_metaData.currentCycle / m_nthTimestep);
    grid->setBlock(rank);
    for (size_t i = 0; i < 3; i++) {
        for (size_t j = 0; j < dims[i]; j++) {
            switch (i) {
            case 0:
                grid->coords(i)[j] = j;
                break;
            case 1:
                grid->coords(i)[j] = m_metaData.currentCycle * dims[i] + j;
                break;
            case 2:
                grid->coords(i)[j] = rank * dims[i] + j;
                break;
            default:
                break;
            }
        }
    }
    CERR << "sending test data for cycle " << m_metaData.currentCycle << endl;
    addObject("AMR_mesh", grid);
}

void insitu::Engine::finalizeInit()     {
    if (!m_moduleInitialized) {
        try {
            getMetaData();
            addPorts();
            getRegisteredGenericCommands();

        } catch (const VistleLibSimExeption& ex) {
            CERR << "finalizeInit failed: " << ex.what() << endl;
            return;
        }
        m_moduleInitialized = true;
    }
}

void insitu::Engine::initializeEngineSocket(const std::string& hostname, int port) {

    boost::system::error_code ec;
    asio::ip::tcp::resolver resolver(m_ioService);
    asio::ip::tcp::resolver::query query(hostname, boost::lexical_cast<std::string>(port));
    m_socket.reset(new boost::asio::ip::tcp::socket(m_ioService));
    asio::ip::tcp::resolver::iterator endpoint_iterator = resolver.resolve(query, ec);
    if (ec) {
        throw EngineExeption("initializeEngineSocket failed to resolve connect socket");
    }
    asio::connect(*m_socket, endpoint_iterator, ec);
    if (ec) {
        throw EngineExeption("initializeEngineSocket failed to connect socket");
    }
}

void insitu::Engine::addObject(const std::string& name, vistle::Object::ptr obj) {
    if (m_sendMessageQueue) {
        DEBUG_CERR << "addObject " << name << " of type " << obj->typeName() << endl;
        vistle::message::AddObject msg(name, obj);
        vistle::message::Buffer buf(msg);

#ifdef MODULE_THREAD
        buf.setSenderId(m_moduleID);
        buf.setRank(m_rank);
#endif
        m_sendMessageQueue->send(buf);
    }


}

Engine::Engine()
{ 

}

Engine::~Engine() {
    m_meshes.clear();
#ifndef MODULE_THREAD
    if (vistle::Shm::isAttached()) {
        vistle::Shm::the().detach();
    }
#endif
    InSituTcpMessage::send(ConnectionClosed{});
    delete m_sendMessageQueue;
    Engine::instance = nullptr;
}







