#include "listenv4v6.h"

#include <boost/asio/ip/v6_only.hpp>

namespace vistle {

bool start_listen(unsigned short port, boost::asio::ip::tcp::acceptor &acceptor_v4, boost::asio::ip::tcp::acceptor &acceptor_v6, boost::system::error_code &ec) {

    using namespace boost::asio;

    boost::system::error_code lec;
    ec = lec;

    acceptor_v4.open(ip::tcp::v4(), lec);
    if (lec == boost::system::errc::address_family_not_supported) {

    } else if (lec) {
        acceptor_v4.close();
        ec = lec;
        return false;
    } else {
        acceptor_v4.bind(ip::tcp::endpoint(ip::tcp::v4(), port), lec);
        if (lec) {
            acceptor_v4.close();
            ec = lec;
            return false;
        }
        acceptor_v4.listen();
    }

    acceptor_v6.open(ip::tcp::v6(), lec);
    if (lec == boost::system::errc::address_family_not_supported) {

    } else if (lec) {
        acceptor_v4.close();
        acceptor_v6.close();
        ec = lec;
        return false;
    } else {
        acceptor_v6.set_option(ip::v6_only(true));
        acceptor_v6.bind(ip::tcp::endpoint(ip::tcp::v6(), port), lec);
        if (lec) {
            acceptor_v4.close();
            acceptor_v6.close();
            ec = lec;
            return false;
        }

        acceptor_v6.listen();
    }

    return true;
}

}
