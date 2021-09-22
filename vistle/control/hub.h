#ifndef VISTLE_HUB_H
#define VISTLE_HUB_H

#include <memory>
#include <atomic>
#include <boost/asio/ip/tcp.hpp>
#include <vistle/core/statetracker.h>
#include <vistle/util/spawnprocess.h>
#include <vistle/util/buffer.h>
#include "uimanager.h"
#include <vistle/net/tunnel.h>
#include <vistle/net/dataproxy.h>
#include <vistle/control/scanmodules.h>

#include "export.h"

namespace vistle {

class V_HUBEXPORT Hub {
public:
    typedef boost::asio::ip::tcp::socket socket;
    typedef boost::asio::ip::tcp::acceptor acceptor;

    static Hub &the();

    Hub(bool inManager = false);
    ~Hub();

    int run();

    bool init(int argc, char *argv[]);
    bool processScript();
    bool dispatch();
    bool sendMessage(std::shared_ptr<socket> sock, const message::Message &msg, const buffer *payload = nullptr);
    unsigned short port() const;
    unsigned short dataPort() const;
    vistle::process_handle launchProcess(const std::vector<std::string> &argv) const;
    const std::string &name() const;

    bool
    handleMessage(const message::Message &msg,
                  std::shared_ptr<boost::asio::ip::tcp::socket> sock = std::shared_ptr<boost::asio::ip::tcp::socket>(),
                  const buffer *payload = nullptr);
    bool sendManager(const message::Message &msg, int hub = message::Id::LocalHub, const buffer *payload = nullptr);
    bool sendMaster(const message::Message &msg, const buffer *payload = nullptr);
    bool sendSlaves(const message::Message &msg, bool returnToSender = false, const buffer *payload = nullptr);
    bool sendHub(int id, const message::Message &msg, const buffer *payload = nullptr);
    bool sendUi(const message::Message &msg, int id = message::Id::Broadcast, const buffer *payload = nullptr);
    bool sendModule(const message::Message &msg, int id, const buffer *payload = nullptr);

    const StateTracker &stateTracker() const;
    StateTracker &stateTracker();

    int idToHub(int id) const;
    int id() const;

private:
    struct Slave;

    message::MessageFactory make;

    bool hubReady();
    bool connectToMaster(const std::string &host, unsigned short port);
    bool startUi(const std::string &uipath, bool replace = false);
    bool startPythonUi();
    bool startServer();
    bool startAccept(std::shared_ptr<acceptor> a);
    void handleWrite(std::shared_ptr<boost::asio::ip::tcp::socket> sock, const boost::system::error_code &error);
    void handleAccept(std::shared_ptr<acceptor> a, std::shared_ptr<boost::asio::ip::tcp::socket> sock,
                      const boost::system::error_code &error);
    void addSocket(std::shared_ptr<boost::asio::ip::tcp::socket> sock,
                   message::Identify::Identity ident = message::Identify::UNKNOWN);
    bool removeSocket(std::shared_ptr<boost::asio::ip::tcp::socket> sock, bool close = true);
    void addClient(std::shared_ptr<boost::asio::ip::tcp::socket> sock);
    void addSlave(const std::string &name, std::shared_ptr<boost::asio::ip::tcp::socket> sock);
    void removeSlave(int id);
    bool removeClient(std::shared_ptr<boost::asio::ip::tcp::socket> sock);
    void slaveReady(Slave &slave);
    bool startCleaner();
    bool processScript(const std::string &filename, bool executeModules);
    bool processStartupScripts();
    void cacheModuleValues(int oldModuleId, int newModuleId);
    void killOldModule(int migratedId);
    void sendInfo(const std::string &s);
    void sendError(const std::string &s);

    bool m_inManager = false;

    unsigned short m_basePort = 31093;
    unsigned short m_port = 0, m_dataPort = 0, m_masterPort = m_basePort;
    std::string m_exposedHost;
    boost::asio::ip::address m_exposedHostAddr;
    std::string m_masterHost;
    boost::asio::io_service m_ioService;
    std::shared_ptr<acceptor> m_acceptorv4, m_acceptorv6;

    std::map<std::shared_ptr<boost::asio::ip::tcp::socket>, message::Identify::Identity> m_sockets;
    std::set<std::shared_ptr<boost::asio::ip::tcp::socket>> m_clients;

    std::shared_ptr<DataProxy> m_dataProxy;
    TunnelManager m_tunnelManager;
    StateTracker m_stateTracker;
    UiManager m_uiManager;
    bool m_hasUi = false;

    std::map<process_handle, int> m_processMap;
    bool m_managerConnected;

    std::string m_prefix;
    std::string m_scriptPath;
    bool m_executeModules = false;
    bool m_quitting = false, m_emergency = false;
    static volatile std::atomic<bool> m_interrupt;
    boost::asio::signal_set m_signals;
    static void signalHandler(const boost::system::error_code &error, int signal_number);

    AvailableMap m_availableModules;
    std::vector<AvailableModule> m_localModules;

    bool m_isMaster;
    std::shared_ptr<boost::asio::ip::tcp::socket> m_masterSocket;
    struct Slave {
        std::shared_ptr<boost::asio::ip::tcp::socket> sock;
        std::string name;
        bool ready = false;
        int id = 0;
    };
    std::map<int, Slave> m_slaves;
    std::vector<Slave *> m_slavesToConnect;
    int m_slaveCount;
    int m_hubId;
    int m_localRanks = -1;
    std::string m_name;
    bool m_ready = false;

    int m_moduleCount;
    message::Type m_traceMessages;

    int m_execCount;

    bool m_barrierActive;
    unsigned m_barrierReached;
    message::uuid_t m_barrierUuid;

    std::string m_statusText;

    std::string sim2FilePath;
    bool m_inSitu = false;

    bool handlePriv(const message::Quit &quit, message::Identify::Identity senderType);
    bool handlePriv(const message::RemoveHub &rm);
    bool handlePriv(const message::Execute &exec);
    bool handlePriv(const message::CancelExecute &cancel);
    bool handlePriv(const message::Barrier &barrier);
    bool handlePriv(const message::BarrierReached &reached);
    bool handlePriv(const message::RequestTunnel &tunnel);
    bool handlePriv(const message::Connect &conn);
    bool handlePriv(const message::Disconnect &disc);
    bool handlePriv(const message::FileQuery &query, const buffer *payload);
    bool handlePriv(const message::FileQueryResult &result, const buffer *payload);
    bool handlePriv(const message::Cover &cover, const buffer *payload);

    bool checkChildProcesses(bool emergency = false);
    bool hasChildProcesses(bool ignoreGui = false);
    void emergencyQuit();
    const AvailableModule *findModule(const AvailableModule::Key &key, int spawnId);
    void spawnModule(const std::string &path, const std::string &name, int spawnId);
    std::vector<message::Buffer> m_queue;
    bool handleQueue();

    void setLoadedFile(const std::string &file);
    void setSessionUrl(const std::string &url);
    void setStatus(const std::string &s, message::UpdateStatus::Importance prio = message::UpdateStatus::Low);
    void clearStatus();

    std::map<int, std::vector<message::Buffer>> m_sendAfterSpawn;

#if BOOST_VERSION >= 106600
    boost::asio::executor_work_guard<boost::asio::io_context::executor_type> m_workGuard;
#else
    std::shared_ptr<boost::asio::io_service::work> m_workGuard;
#endif
    std::thread m_ioThread;
};

} // namespace vistle
#endif
