#include "random-cache-strategy.hpp"

#include <iostream>

#include "ns3/core-module.h"
#include "ns3/ndnSIM/NFD/daemon/fw/algorithm.hpp"
#include "ns3/scheduler.h"
#include "ns3/simulator.h"

namespace nfd {
namespace fw {

NFD_LOG_INIT(RandomCacheStrategy);
NFD_REGISTER_STRATEGY(RandomCacheStrategy);

RandomCacheStrategy::RandomCacheStrategy(Forwarder& forwarder, const Name& name)
    : Strategy(forwarder), fwd(forwarder), boundary2app(0) {
    ParsedInstanceName parsed = parseInstanceName(name);
    std::cout << "strategy ptr = " << this << '\n';
    if (parsed.parameters.empty()) {
        std::cout << "\tno parameter\n";
        std::cout << "\tboundary2app = " << boundary2app << std::endl;
    } else {
        auto para = std::stod(parsed.parameters[0].toUri());
        std::cout << "\tparameter is " << para << '\n';
        boundary2app = MAX64 * para;
        std::cout << "\tboundary2app = " << boundary2app << std::endl;
    }
    if (parsed.version &&
        *parsed.version != getStrategyName()[-1].toVersion()) {
        NDN_THROW(std::invalid_argument(
            "RandomCacheStrategy does not support version " +
            to_string(*parsed.version)));
    }
    this->setInstanceName(makeInstanceName(name, getStrategyName()));

    //::ns3::Simulator::Schedule(::ns3::MilliSeconds(500), ::ns3::MakeCallback());
}

const Name& RandomCacheStrategy::getStrategyName() {
    static Name strategyName(
        "/localhost/nfd/strategy/random-cache-strategy/%FD%01");
    return strategyName;
}

void RandomCacheStrategy::afterReceiveInterest(
    const FaceEndpoint& ingress, const Interest& interest,
                                            const shared_ptr<pit::Entry>& pitEntry){
        
    if (hasPendingOutRecords(*pitEntry)) {
        // not a new Interest, don't forward
        return;
    }

    auto regular = [&, this](const fib::NextHopList& nexthops) {
        // std::cout << "regular\n";
        // for (const auto& nexthop : nexthops) {
        //     std::cout << this << ' ' << nexthop.getFace().getRemoteUri()
        //               << '\n';
        // }
        // std::cout << std::endl;

        for (const auto& nexthop : nexthops) {
            Face& outFace = nexthop.getFace();
            if (isNextHopEligible(ingress.face, interest, nexthop, pitEntry,
                                  false, time::steady_clock::time_point())) {
                this->sendInterest(pitEntry, FaceEndpoint(outFace, 0), interest);
                return;
            }
        }
        this->rejectPendingInterest(pitEntry);
    };

    auto random = [&, this](const fib::NextHopList& nexthops) {
        // std::cout << "random\n";
        // for (const auto& nexthop : nexthops) {
        //     std::cout << this << ' ' << nexthop.getFace().getRemoteUri()
        //               << '\n';
        // }
        fib::NextHopList canforwardhops;
        for (const auto& nexthop : nexthops) {
            Face& outFace = nexthop.getFace();
            if (isNextHopEligible(ingress.face, interest, nexthop, pitEntry,
                                  false, time::steady_clock::time_point())) {
                canforwardhops.emplace_back(nexthop);
            }
        }
        // std::cout << "canforward\n";
        // for (const auto& nexthop : canforwardhops) {
        //     std::cout << this << ' ' << nexthop.getFace().getRemoteUri()
        //               << '\n';
        // }
        // std::cout << std::endl;

        if (!canforwardhops.empty()) {
            std::random_device rdev;
            std::mt19937 reng(rdev());
            std::uniform_int_distribution<int> unir(0,
                                                    canforwardhops.size() - 1);
            this->sendInterest(pitEntry, FaceEndpoint(canforwardhops[unir(reng)].getFace(),0),
                               interest);
            return;
        }

        this->rejectPendingInterest(pitEntry);
    };

    const fib::Entry& fibEntry = this->lookupFib(*pitEntry);
    const fib::NextHopList& nexthops_origin = fibEntry.getNextHops();
    if (generateWord64() >= boundary2app) {
        fib::NextHopList nexthops_manip;
        for (auto& it : nexthops_origin) {
            if (it.getFace().getRemoteUri().getScheme() != "appFace") {
                nexthops_manip.emplace_back(it);
            }
        }
        random(nexthops_manip);
    } else {
        regular(nexthops_origin);
    }
}

// void RandomCacheStrategy::runevery1s() {
//     boundary2app = MAX64 * prob;

//     if()
//     ::ns3::Simulator::Schedule(::ns3::MilliSeconds(500), ::ns3::MakeCallback(&RandomCacheStrategy::runevery1s,this))
// }

}  // namespace fw
}  // namespace nfd