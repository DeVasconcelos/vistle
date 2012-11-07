/*
 * Visualization Testing Laboratory for Exascale Computing (VISTLE)
 */
#ifndef _WIN32
#include <sys/socket.h>
#include <sys/select.h>
#include <netinet/in.h>
#include <unistd.h>
#else
#include<Winsock2.h>
#pragma comment(lib, "Ws2_32.lib")
#define usleep Sleep
#define close closesocket
typedef int socklen_t;
#endif

#include <mpi.h>

#include <sys/types.h>

#include <stdlib.h>
#include <stdio.h>

#include <sstream>
#include <iostream>
#include <iomanip>

#include "message.h"
#include "messagequeue.h"
#include "object.h"
#include "shm.h"

#include "communicator.h"

//#define CLIENTSOCKET

using namespace boost::interprocess;

namespace vistle {


#ifdef CLIENTSOCKET
static int acceptClient();
#endif

static void splitpath(const std::string &value, std::vector<std::string> *components)
{
#ifdef _WIN32
   const char *sep = ";";
#else
   const char *sep = ":";
#endif

   std::string::size_type begin = 0;
   do {
      std::string::size_type end = value.find(sep, begin);

      std::string c;
      if (end != std::string::npos) {
         c = value.substr(begin, end-begin);
         ++end;
      } else {
         c = value.substr(begin);
      }
      begin = end;

      if (!c.empty())
         components->push_back(c);
   } while(begin != std::string::npos);
}

static std::string getbindir(int argc, char *argv[]) {

   char *wd = getcwd(NULL, 0);
   if (!wd) {

      std::cerr << "failed to determine working directory: " << strerror(errno) << std::endl;
      exit(1);
   }
   std::string cwd(wd);
   free(wd);

   // determine complete path to executable
   std::string executable;
#ifdef _WIN32
   char buf[2000];
   DWORD sz = GetModuleFileName(NULL, buf, sizeof(buf));
   if (sz != 0) {

      executable = buf;
   } else {

      std::cerr << "getbindir(): GetModuleFileName failed - error: " << GetLastError() << std::endl;
   }
#else
   char buf[PATH_MAX];
   ssize_t len = readlink("/proc/self/exe", buf, sizeof(buf));
   if(len != -1) {

      executable = std::string(buf, len);
   } else if (argc >= 1) {

      bool found = false;
      if (!strchr(argv[0], '/')) {

         if(const char *path = getenv("PATH")) {
            std::vector<std::string> components;
            splitpath(path, &components);
            for(size_t i=0; i<components.size(); ++i) {

               std::string component = components[i];
               if (component[0] != '/')
                  component = cwd + "/" + component;

               DIR *dir = opendir(component.c_str());
               if (!dir) {
                  std::cerr << "failed to open directory " << component << ": " << strerror(errno) << std::endl;
                  continue;
               }

               while (struct dirent *ent = readdir(dir)) {
                  if (!strcmp(ent->d_name, argv[0])) {
                     found = true;
                     break;
                  }
               }

               if (found) {
                  executable = component + "/" + argv[0];
                  break;
               }
            }
         }
      } else if (argv[0][0] != '/') {
         executable = cwd + "/" + argv[0];
      } else {
         executable = argv[0];
      }
   }
#endif

   // guess vistle bin directory
   if (!executable.empty()) {

      std::string dir = executable;
#ifdef _WIN32
      std::string::size_type idx = dir.find_last_of("\\/");
#else
      std::string::size_type idx = dir.rfind('/');
#endif
      if (idx == std::string::npos) {

         dir = cwd;
      } else {

         dir = executable.substr(0, idx);
      }

#ifdef DEBUG
      std::cerr << "vistle bin directory determined to be " << dir << std::endl;
#endif
      return dir;
   }

   return std::string();
}

Communicator::Communicator(int argc, char *argv[], int r, int s)
   : rank(r), size(s), socketBuffer(NULL), clientSocket(-1), moduleID(0),
     mpiReceiveBuffer(NULL), mpiMessageSize(0) {

   m_bindir = getbindir(argc, argv);

   socketBuffer = new unsigned char[64];
   mpiReceiveBuffer = new char[message::Message::MESSAGE_SIZE];

   // post requests for length of next MPI message
   MPI_Irecv(&mpiMessageSize, 1, MPI_INT, MPI_ANY_SOURCE, 0,
             MPI_COMM_WORLD, &request);
   MPI_Barrier(MPI_COMM_WORLD);

#ifdef CLIENTSOCKET
   if (rank == 0)
      clientSocket = acceptClient();
#endif
}

int Communicator::getRank() const {

   return rank;
}

int Communicator::getSize() const {

   return size;
}

bool Communicator::dispatch() {

   bool done = false;

   if (rank == 0 && clientSocket != -1) {
      // poll socket
      fd_set set;
      FD_ZERO(&set);
      FD_SET(clientSocket, &set);

      struct timeval t = { 0, 0 };

      select(clientSocket + 1, &set, NULL, NULL, &t);

      if (FD_ISSET(clientSocket, &set)) {

         message::Message *message = NULL;

         int r = recv(clientSocket, (char *) socketBuffer, 1,0);
         if (r == 1) {

            if (socketBuffer[0] == 'q')
               message = new message::Quit(0, rank);

            else if (socketBuffer[0] != '\r' && socketBuffer[0] != '\n')
               message = new message::Debug(0, rank, socketBuffer[0]);

            // Broadcast message to other MPI partitions
            //   - handle message, delete message
            if (message) {

               MPI_Request s;
               for (int index = 0; index < size; index ++)
                  if (index != rank)
                     MPI_Isend(&(message->size), 1, MPI_INT, index, 0,
                               MPI_COMM_WORLD, &s);

               MPI_Bcast(message, message->size, MPI_BYTE, 0, MPI_COMM_WORLD);

               if (!handleMessage(*message))
                  done = true;

               delete message;
            }
         }
      }
   }

   int flag;
   MPI_Status status;

   // test for message size from another MPI node
   //    - receive actual message from broadcast
   //    - handle message
   //    - post another MPI receive for size of next message
   MPI_Test(&request, &flag, &status);

   if (flag) {
      if (status.MPI_TAG == 0) {

         MPI_Bcast(mpiReceiveBuffer, mpiMessageSize, MPI_BYTE,
                   status.MPI_SOURCE, MPI_COMM_WORLD);

         message::Message *message = (message::Message *) mpiReceiveBuffer;
#if 0
         printf("[%02d] message from [%02d] message type %d size %d\n",
                rank, status.MPI_SOURCE, message->getType(), mpiMessageSize);
#endif
         if (!handleMessage(*message))
            done = true;

         MPI_Irecv(&mpiMessageSize, 1, MPI_INT, MPI_ANY_SOURCE, 0,
                   MPI_COMM_WORLD, &request);
      }
   }

   // test for messages from modules
   size_t msgSize;
   unsigned int priority;
   char msgRecvBuf[vistle::message::Message::MESSAGE_SIZE];

   std::map<int, message::MessageQueue *>::iterator i;
   for (i = receiveMessageQueue.begin(); i != receiveMessageQueue.end(); ){

      bool moduleExit = false;
      try {
         bool received =
            i->second->getMessageQueue().try_receive(
                                        (void *) msgRecvBuf,
                                        vistle::message::Message::MESSAGE_SIZE,
                                        msgSize, priority);

         if (received) {
            moduleExit = !handleMessage(*(message::Message *) msgRecvBuf);

            if (moduleExit) {

               std::map<int, message::MessageQueue *>::iterator si =
                  sendMessageQueue.find(i->first);
               if (si != sendMessageQueue.end()) {
                  delete si->second;
                  sendMessageQueue.erase(si);
               }
               receiveMessageQueue.erase(i++);
            }
         }
      } catch (interprocess_exception &ex) {
         std::cerr << "comm [" << rank << "/" << size << "] receive mq "
                   << ex.what() << std::endl;
         exit(-1);
      }
      if (!moduleExit)
         i++;
   }

   return done;
}


bool Communicator::handleMessage(const message::Message &message) {

   switch (message.getType()) {

      case message::Message::DEBUG: {

         const message::Debug &debug =
            static_cast<const message::Debug &>(message);
         std::cout << "comm [" << rank << "/" << size << "] Debug ["
                   << debug.getCharacter() << "]" << std::endl;
         break;
      }

      case message::Message::QUIT: {

         const message::Quit &quit =
            static_cast<const message::Quit &>(message);
         (void) quit;
         return false;
         break;
      }

      case message::Message::SPAWN: {

         const message::Spawn &spawn =
            static_cast<const message::Spawn &>(message);
         int moduleID = spawn.getSpawnID();

         std::string name = m_bindir + "/" + spawn.getName();

         std::stringstream modID;
         modID << moduleID;

         std::string smqName =
            message::MessageQueue::createName("smq", moduleID, rank);
         std::string rmqName =
            message::MessageQueue::createName("rmq", moduleID, rank);

         try {
            sendMessageQueue[moduleID] =
               message::MessageQueue::create(smqName);
            receiveMessageQueue[moduleID] =
               message::MessageQueue::create(rmqName);
         } catch (interprocess_exception &ex) {

            std::cerr << "comm [" << rank << "/" << size << "] spawn mq "
                      << ex.what() << std::endl;
            exit(-1);
         }

         MPI_Comm interComm;
         char *argv[3] = { strdup(Shm::the().name().c_str()), strdup(modID.str().c_str()), NULL };

         MPI_Comm_spawn(strdup(name.c_str()), argv, size, MPI_INFO_NULL,
                        0, MPI_COMM_WORLD, &interComm, MPI_ERRCODES_IGNORE);

         break;
      }

      case message::Message::CONNECT: {

         const message::Connect &connect =
            static_cast<const message::Connect &>(message);
         portManager.addConnection(connect.getModuleA(),
                                   connect.getPortAName(),
                                   connect.getModuleB(),
                                   connect.getPortBName());
         break;
      }

      case message::Message::NEWOBJECT: {

         /*
         const message::NewObject *newObject =
            static_cast<const message::NewObject *>(message);
         vistle::Object *object = (vistle::Object *)
            vistle::Shm::the().shm().get_address_from_handle(newObject.getHandle());

         std::cout << "comm [" << rank << "/" << size << "] NewObject ["
                   << newObject.getHandle() << "] type ["
                   << object.getType() << "] from module ["
                   << newObject.getModuleID() << "]" << std::endl;
         */
         break;
      }

      case message::Message::MODULEEXIT: {

         const message::ModuleExit &moduleExit =
            static_cast<const message::ModuleExit &>(message);
         int mod = moduleExit.getModuleID();

         std::cout << "comm [" << rank << "/" << size << "] Module ["
                   << mod << "] quit" << std::endl;

         return false;
         break;
      }

      case message::Message::COMPUTE: {

         const message::Compute &comp =
            static_cast<const message::Compute &>(message);
         std::map<int, message::MessageQueue *>::iterator i
            = sendMessageQueue.find(comp.getModule());
         if (i != sendMessageQueue.end())
            i->second->getMessageQueue().send(&comp, sizeof(comp), 0);
         break;
      }

      case message::Message::CREATEINPUTPORT: {

         const message::CreateInputPort &m =
            static_cast<const message::CreateInputPort &>(message);
         portManager.addPort(m.getModuleID(), m.getName(),
                             Port::INPUT);
         break;
      }

      case message::Message::CREATEOUTPUTPORT: {

         const message::CreateOutputPort &m =
            static_cast<const message::CreateOutputPort &>(message);
         portManager.addPort(m.getModuleID(), m.getName(),
                             Port::OUTPUT);
         break;
      }

      case message::Message::ADDOBJECT: {

         const message::AddObject &m =
            static_cast<const message::AddObject &>(message);
         Object::const_ptr obj = m.takeObject();
#if 0
         std::cout << "Module " << m.getModuleID() << ": "
                   << "AddObject " << m.getHandle() << " (" << obj->getName() << ")"
                   << " ref " << obj->refcount()
                   << " to port " << m.getPortName() << std::endl;
#endif

         Port *port = portManager.getPort(m.getModuleID(),
                                          m.getPortName());
         if (port) {
            const std::vector<const Port *> *list =
               portManager.getConnectionList(port);

            std::vector<const Port *>::const_iterator pi;
            for (pi = list->begin(); pi != list->end(); pi ++) {

               std::map<int, message::MessageQueue *>::iterator mi =
                  sendMessageQueue.find((*pi)->getModuleID());
               if (mi != sendMessageQueue.end()) {
                  const message::AddObject a(m.getModuleID(), m.getRank(),
                                             (*pi)->getName(), obj);
                  const message::Compute c(moduleID, rank,
                                           (*pi)->getModuleID());

                  mi->second->getMessageQueue().send(&a, sizeof(a), 0);
                  mi->second->getMessageQueue().send(&c, sizeof(c), 0);
               }
            }
         }
         else
            std::cout << "comm [" << rank << "/" << size << "] Addbject ["
                      << m.getHandle() << "] to port ["
                      << m.getPortName() << "]: port not found" << std::endl;

         break;
      }

      case message::Message::SETFILEPARAMETER: {

         const message::SetFileParameter &m =
            static_cast<const message::SetFileParameter &>(message);

         if (m.getModuleID() != m.getModule()) {
            // message to module
            std::map<int, message::MessageQueue *>::iterator i
               = sendMessageQueue.find(m.getModule());
            if (i != sendMessageQueue.end())
               i->second->getMessageQueue().send(&m, m.getSize(), 0);
         }
         break;
      }

      case message::Message::SETFLOATPARAMETER: {

         const message::SetFloatParameter &m =
            static_cast<const message::SetFloatParameter &>(message);

         if (m.getModuleID() != m.getModule()) {
            // message to module
            std::map<int, message::MessageQueue *>::iterator i
               = sendMessageQueue.find(m.getModule());
            if (i != sendMessageQueue.end())
               i->second->getMessageQueue().send(&m, m.getSize(), 0);
         }
         break;
      }

      case message::Message::SETINTPARAMETER: {

         const message::SetIntParameter &m =
            static_cast<const message::SetIntParameter &>(message);

         if (m.getModuleID() != m.getModule()) {
            // message to module
            std::map<int, message::MessageQueue *>::iterator i
               = sendMessageQueue.find(m.getModule());
            if (i != sendMessageQueue.end())
               i->second->getMessageQueue().send(&m, m.getSize(), 0);
         }
         break;

      }

      case message::Message::SETVECTORPARAMETER: {

         const message::SetVectorParameter &m =
            static_cast<const message::SetVectorParameter &>(message);

         if (m.getModuleID() != m.getModule()) {
            // message to module
            std::map<int, message::MessageQueue *>::iterator i
               = sendMessageQueue.find(m.getModule());
            if (i != sendMessageQueue.end())
               i->second->getMessageQueue().send(&m, m.getSize(), 0);
         }
         break;
      }

      default:
         break;

   }

   return true;
}

Communicator::~Communicator() {

   message::Quit quit(0, rank);
   std::map<int, message::MessageQueue *>::iterator i;

   for (i = sendMessageQueue.begin(); i != sendMessageQueue.end(); i ++)
      i->second->getMessageQueue().send(&quit, sizeof(quit), 1);

   // receive all ModuleExit messages from modules
   // retry for some time, modules that don't answer might have crashed
   int retries = 10000;
   while (sendMessageQueue.size() > 0 && --retries >= 0) {
      dispatch();
      usleep(1000);
   }

   if (clientSocket != -1)
      close(clientSocket);

   if (size > 1) {
      int dummy;
      MPI_Request s;
      MPI_Isend(&dummy, 1, MPI_INT, (rank + 1) % size, 0, MPI_COMM_WORLD, &s);
      MPI_Wait(&request, MPI_STATUS_IGNORE);
      MPI_Wait(&s, MPI_STATUS_IGNORE);
   }
   MPI_Barrier(MPI_COMM_WORLD);
}

#ifdef CLIENTSOCKET
int acceptClient() {

   int server = socket(AF_INET, SOCK_STREAM, 0);
   int reuse = 1;

   setsockopt(server, SOL_SOCKET, SO_REUSEADDR, (char *) &reuse,
              sizeof(reuse));

   struct sockaddr_in serv, addr;
   serv.sin_family = AF_INET;
   serv.sin_addr.s_addr = INADDR_ANY;
   serv.sin_port = htons(8192);

   if (bind(server, (struct sockaddr *) &serv, sizeof(serv)) < 0) {
      perror("bind error");
      exit(1);
   }
   listen(server, 0);

   socklen_t len = sizeof(addr);
   int client = accept(server, (struct sockaddr *) &addr, &len);
   close(server);

   return client;
}
#endif

} // namespace vistle
