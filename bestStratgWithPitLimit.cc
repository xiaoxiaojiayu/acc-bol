/* -*- Mode:C++; c-file-style:"gnu"; indent-tabs-mode:nil; -*- */
/*
 * Copyright (c) 2014-2019,  Regents of the University of California,
 *                           Arizona Board of Regents,
 *                           Colorado State University,
 *                           University Pierre & Marie Curie, Sorbonne University,
 *                           Washington University in St. Louis,
 *                           Beijing Institute of Technology,
 *                           The University of Memphis.
 *
 * This file is part of NFD (Named Data Networking Forwarding Daemon).
 * See AUTHORS.md for complete list of NFD authors and contributors.
 *
 * NFD is free software: you can redistribute it and/or modify it under the terms
 * of the GNU General Public License as published by the Free Software Foundation,
 * either version 3 of the License, or (at your option) any later version.
 *
 * NFD is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 * PURPOSE.  See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * NFD, e.g., in COPYING.md file.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "bestStratgWithPitLimit.h"
#include "ns3/ndnSIM/NFD/daemon/fw/algorithm.hpp"
#include "ns3/ndnSIM/NFD/daemon/common/logger.hpp"
#include "ns3/ndnSIM/NFD/core/common.hpp"

namespace nfd {
    namespace fw {

        NFD_LOG_INIT(BestRouteStrategy2WithPitLimit);
        NFD_REGISTER_STRATEGY(BestRouteStrategy2WithPitLimit);

        const time::milliseconds BestRouteStrategy2WithPitLimit::RETX_SUPPRESSION_INITIAL(10);
        const time::milliseconds BestRouteStrategy2WithPitLimit::RETX_SUPPRESSION_MAX(250);

        BestRouteStrategy2WithPitLimit::BestRouteStrategy2WithPitLimit(Forwarder& forwarder, const Name& name)
            : Strategy(forwarder)
            , ProcessNackTraits(this)
            , m_retxSuppression(RETX_SUPPRESSION_INITIAL, RetxSuppressionExponential::DEFAULT_MULTIPLIER, RETX_SUPPRESSION_MAX)
            , pitMaxSize(230), curPitSize(0), firstrun(true)
            , forwarderAcc(forwarder)
        {
            ParsedInstanceName parsed = parseInstanceName(name);
            if (!parsed.parameters.empty()) {
                NDN_THROW(std::invalid_argument("BestRouteStrategy2WithConges does not accept parameters"));
            }
            if (parsed.version && *parsed.version != getStrategyName()[-1].toVersion()) {
                NDN_THROW(std::invalid_argument("BestRouteStrategy2WithConges does not support version "
                                                + to_string(*parsed.version)));
            }
            this->setInstanceName(makeInstanceName(name, getStrategyName()));
        }

        const Name& BestRouteStrategy2WithPitLimit::getStrategyName()
        {
            static Name strategyName("/localhost/nfd/strategy/best-route2-conges/%FD%01");
            return strategyName;
        }



        void BestRouteStrategy2WithPitLimit::afterReceiveInterest(const FaceEndpoint& ingress, const Interest& interest,
                                                                  const shared_ptr<pit::Entry>& pitEntry)
        {
            ns3::Simulator::Now().GetSeconds();
            if (firstrun) {
                curPitSize = forwarderAcc.getPit().size();
                firstrun = false;
            }
            curPitSize = forwarderAcc.getPit().size();//0.75 * curPitSize + 0.25 * forwarderAcc.getPit().size();
            //std::cout << "curPitSize: " << curPitSize << std::endl;
            if (curPitSize > pitMaxSize) {
                lp::NackHeader nackHeader;
                nackHeader.setReason(lp::NackReason::CONGESTION);
                this->sendNack(pitEntry, ingress, nackHeader);
                this->rejectPendingInterest(pitEntry);
                return;
            }
            else if (curPitSize > 0.5 * pitMaxSize) {
                lp::NackHeader nackHeader;
                nackHeader.setReason(lp::NackReason::BUSY);
                this->sendNack(pitEntry, ingress, nackHeader);
            }
            else {
                lp::NackHeader nackHeader;
                nackHeader.setReason(lp::NackReason::FREE);
                this->sendNack(pitEntry, ingress, nackHeader);
            }

            RetxSuppressionResult suppression = m_retxSuppression.decidePerPitEntry(*pitEntry);
            if (suppression == RetxSuppressionResult::SUPPRESS) {
                NFD_LOG_DEBUG(interest << " from=" << ingress << " suppressed");
                return;
            }

            const fib::Entry& fibEntry = this->lookupFib(*pitEntry);
            const fib::NextHopList& nexthops = fibEntry.getNextHops();
            auto it = nexthops.end();

            if (suppression == RetxSuppressionResult::NEW) {
                // forward to nexthop with lowest cost except downstream
                it = std::find_if(nexthops.begin(), nexthops.end(), [&](const auto& nexthop) {
                    return isNextHopEligible(ingress.face, interest, nexthop, pitEntry);
                                  });

                if (it == nexthops.end()) {
                    NFD_LOG_DEBUG(interest << " from=" << ingress << " noNextHop");

                    lp::NackHeader nackHeader;
                    nackHeader.setReason(lp::NackReason::NO_ROUTE);
                    this->sendNack(pitEntry, ingress, nackHeader);

                    this->rejectPendingInterest(pitEntry);
                    return;
                }

                auto egress = FaceEndpoint(it->getFace(), 0);
                //if (egress.face.getTransport()->getSendQueueLength() >= 0)
                //    std::cout << "cur queue lenth of out face: " << egress.face.getTransport()->getSendQueueLength()
                //    << "/" << egress.face.getTransport()->getSendQueueCapacity() << std::endl;

                NFD_LOG_DEBUG(interest << " from=" << ingress << " newPitEntry-to=" << egress);
                this->sendInterest(pitEntry, egress, interest);
                return;
            }

            // find an unused upstream with lowest cost except downstream
            it = std::find_if(nexthops.begin(), nexthops.end(), [&](const auto& nexthop) {
                return isNextHopEligible(ingress.face, interest, nexthop, pitEntry, true, time::steady_clock::now());
                              });

            if (it != nexthops.end()) {
                auto egress = FaceEndpoint(it->getFace(), 0);
                //if (egress.face.getTransport()->getSendQueueLength() >= 0)
                //    std::cout << "cur queue lenth of out face: " << egress.face.getTransport()->getSendQueueLength()
                //    << "/" << egress.face.getTransport()->getSendQueueCapacity() << std::endl;

                this->sendInterest(pitEntry, egress, interest);
                NFD_LOG_DEBUG(interest << " from=" << ingress << " retransmit-unused-to=" << egress);
                return;
            }

            // find an eligible upstream that is used earliest
            it = findEligibleNextHopWithEarliestOutRecord(ingress.face, interest, nexthops, pitEntry);
            if (it == nexthops.end()) {
                NFD_LOG_DEBUG(interest << " from=" << ingress << " retransmitNoNextHop");
            }
            else {
                auto egress = FaceEndpoint(it->getFace(), 0);
                //if (egress.face.getTransport()->getSendQueueLength() >= 0)
                //    std::cout << "cur queue lenth of out face: " << egress.face.getTransport()->getSendQueueLength()
                //    << "/" << egress.face.getTransport()->getSendQueueCapacity() << std::endl;
                this->sendInterest(pitEntry, egress, interest);
                NFD_LOG_DEBUG(interest << " from=" << ingress << " retransmit-retry-to=" << egress);
            }
        }

        void BestRouteStrategy2WithPitLimit::afterReceiveNack(const FaceEndpoint& ingress, const lp::Nack& nack,
                                                              const shared_ptr<pit::Entry>& pitEntry)
        {
            this->processNack(ingress.face, nack, pitEntry);
        }

    } // namespace fw
} // namespace nfd
