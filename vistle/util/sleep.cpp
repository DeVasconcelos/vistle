#include "sleep.h"
#include "sysdep.h"
#include <map>
#include <iostream>
#include <mutex>

#ifndef _WIN32
#include <unistd.h>
#include <sched.h>
#endif

namespace vistle {

bool adaptive_wait(bool work, const void *client)
{
    static std::mutex protect;
    static std::map<const void *, long> idleMap;
    const long Sec = 1000000; // 1 s
    const long MinDelay = Sec / 10000;
    const long MaxDelay = Sec / 100;

    protect.lock();
    auto it = idleMap.find(client);
    if (it == idleMap.end())
        it = idleMap.insert(std::make_pair(client, 0)).first;

    auto &idle = it->second;
    protect.unlock();

    if (work) {
        idle = 0;
        return false;
    }


    long delay = (long)((float(idle) / Sec) * Sec * 0.1f);
    if (delay < MinDelay)
        delay = MinDelay;
    if (delay > MaxDelay)
        delay = MaxDelay;

    idle += delay;

    if (idle > delay) {
#if 0
      if (delay == MinDelay)
         std::cerr << "." << std::flush;
      else if (delay == MaxDelay)
         std::cerr << "O" << std::flush;
      else
         std::cerr << "o" << std::flush;
#endif
        if (delay < Sec) {
            //std::cerr << "usleep " << delay << std::endl;
            usleep(delay);
        } else {
            sleep(delay / Sec);
            //std::cerr << "sleep " << delay/Sec << std::endl;
        }
    } else {
#ifdef _POSIX_PRIORITY_SCHEDULING
        sched_yield();
#endif
    }

    return true;
}

} // namespace vistle
