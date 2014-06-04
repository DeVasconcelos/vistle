#include <sstream>
#include <iomanip>

#include "message.h"
#include "messagequeue.h"
#include "shm.h"

using namespace boost::interprocess;

namespace vistle {
namespace message {

std::string MessageQueue::createName(const char * prefix,
                                     const int moduleID, const int rank) {

   std::stringstream mqID;
   mqID << "mq_" << Shm::the().name() << "_" << prefix << "_" << std::setw(8) << std::setfill('0') << moduleID
        << "_" << std::setw(8) << std::setfill('0') << rank;

   return mqID.str();
}

MessageQueue * MessageQueue::create(const std::string & n) {

   {
      std::ofstream f;
      f.open(Shm::shmIdFilename().c_str(), std::ios::app);
      f << n << std::endl;
   }

   message_queue::remove(n.c_str());
   return new MessageQueue(n, create_only);
}

MessageQueue * MessageQueue::open(const std::string & n) {

   auto ret = new MessageQueue(n, open_only);
   message_queue::remove(n.c_str());
   return ret;
}

MessageQueue::MessageQueue(const std::string & n, create_only_t)
   : m_name(n),  m_mq(create_only, m_name.c_str(), 256 /* num msg */,
                  message::Message::MESSAGE_SIZE)
{
}

MessageQueue::MessageQueue(const std::string & n, open_only_t)
   : m_name(n), m_mq(open_only, m_name.c_str())
{
}

MessageQueue::~MessageQueue() {

   message_queue::remove(m_name.c_str());
}

const std::string & MessageQueue::getName() const {

   return m_name;
}

void MessageQueue::send(const Message &msg) {

   Buffer buf(msg);
   m_mq.send(&buf.msg, message::Message::MESSAGE_SIZE, 0);
}

void MessageQueue::receive(Message &msg) {

   size_t recvSize = 0;
   unsigned priority = 0;
   m_mq.receive(&msg, message::Message::MESSAGE_SIZE, recvSize, priority);
   vassert(recvSize == message::Message::MESSAGE_SIZE);
}

bool MessageQueue::tryReceive(Message &msg) {

   size_t recvSize = 0;
   unsigned priority = 0;
   bool result = m_mq.try_receive(&msg, message::Message::MESSAGE_SIZE, recvSize, priority);
   if (result) {
      vassert(recvSize == message::Message::MESSAGE_SIZE);
   }
   return result;
}

size_t MessageQueue::getNumMessages() {

   return m_mq.get_num_msg();
}

} // namespace message
} // namespace vistle
