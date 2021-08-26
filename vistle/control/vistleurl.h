#ifndef VISTLE_URL_H
#define VISTLE_URL_H

#include "export.h"
#include <string>

namespace vistle {

class ConnectionData {
public:
    bool master = false;
    std::string host;
    unsigned short port = 0;
    std::string hex_key;
    std::string kind;
};


class VistleUrl {

public:
    static bool parse(std::string url, ConnectionData &data);
    static std::string create(std::string host, unsigned short port, std::string key);
    static std::string create(const ConnectionData &data);
};

}
#endif //VISTLE_URL_H
