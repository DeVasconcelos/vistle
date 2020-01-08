#ifndef VISTLE_TUNNEL_H
#define VISTLE_TUNNEL_H

#include <memory>
#include <boost/asio.hpp>
#include <thread>

#include <core/message.h>
#include <core/messages.h>

#include "export.h"

namespace vistle {

class TunnelManager;
class TunnelStream;

//! one tunnel waiting for connections
class V_NETEXPORT Tunnel {

   typedef boost::asio::ip::tcp::socket socket;
   typedef boost::asio::ip::tcp::acceptor acceptor;
   typedef boost::asio::ip::address address;
public:
   Tunnel(TunnelManager &tunnel, unsigned short listenPort,
         address destAddr, unsigned short destPort);
   ~Tunnel();
   void startAccept(std::shared_ptr<Tunnel> self);
   void shutdown();

   void cleanUp();
   const address &destAddr() const;
   unsigned short destPort() const;

private:
   TunnelManager &m_manager;
   address m_destAddr;
   unsigned short m_destPort;

   acceptor m_acceptorv4, m_acceptorv6;
   std::map<acceptor *, std::shared_ptr<socket>> m_listeningSocket;
   std::vector<std::weak_ptr<TunnelStream>> m_streams;
   void startAccept(acceptor &a, std::shared_ptr<Tunnel> self);
   void handleAccept(acceptor &a, std::shared_ptr<Tunnel> self, const boost::system::error_code &error);
   void handleConnect(std::shared_ptr<Tunnel> self, std::shared_ptr<boost::asio::ip::tcp::socket> sock0, std::shared_ptr<boost::asio::ip::tcp::socket> sock1, const boost::system::error_code &error);
};

//! a single established connection being tunneled
class V_NETEXPORT TunnelStream {

   typedef boost::asio::ip::tcp::socket socket;
 public:
   TunnelStream(std::shared_ptr<boost::asio::ip::tcp::socket> sock0, std::shared_ptr<boost::asio::ip::tcp::socket> sock1);
   ~TunnelStream();
   void start(std::shared_ptr<TunnelStream> self);
   void destroy();
   void close();

 private:
   std::vector<std::vector<char>> m_buf;
   std::vector<std::shared_ptr<socket>> m_sock;
   void handleRead(std::shared_ptr<TunnelStream> self, size_t sockIdx, boost::system::error_code ec, size_t length);
   void handleWrite(std::shared_ptr<TunnelStream> self, size_t sockIdx, boost::system::error_code ec);
};

//! manage tunnel creation and destruction
class V_NETEXPORT TunnelManager {

 public:
   typedef boost::asio::ip::tcp::socket socket;
   typedef boost::asio::ip::tcp::acceptor acceptor;
   typedef boost::asio::io_service io_service;

   TunnelManager();
   ~TunnelManager();
   io_service &io();

   bool processRequest(const message::RequestTunnel &msg);
   void cleanUp();

 private:
   bool addTunnel(const message::RequestTunnel &msg);
   bool removeTunnel(const message::RequestTunnel &msg);
   void startThread();
   io_service m_io;
   std::map<unsigned short, std::shared_ptr<Tunnel>> m_tunnels;
   std::vector<std::thread> m_threads;
};

}
#endif
   

