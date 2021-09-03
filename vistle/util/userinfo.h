#ifndef VISTLE_USERINFO_H
#define VISTLE_USERINFO_H

#include "export.h"
#include <string>

namespace vistle {

V_UTILEXPORT std::string getLoginName();
V_UTILEXPORT std::string getRealName();
V_UTILEXPORT unsigned long getUserId();

}
#endif
