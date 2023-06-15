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

#include "bestStratgWithECP.h"
#include "ns3/ndnSIM/NFD/daemon/fw/algorithm.hpp"
#include "ns3/ndnSIM/NFD/daemon/common/logger.hpp"
#include "ns3/ndnSIM/NFD/core/common.hpp"

namespace nfd {
    namespace fw {

        NFD_LOG_INIT(BestRouteStrategy2WithECP);
        NFD_REGISTER_STRATEGY(BestRouteStrategy2WithECP);

        const time::milliseconds BestRouteStrategy2WithECP::RETX_SUPPRESSION_INITIAL(10);
        const time::milliseconds BestRouteStrategy2WithECP::RETX_SUPPRESSION_MAX(250);

        BestRouteStrategy2WithECP::BestRouteStrategy2WithECP(Forwarder& forwarder, const Name& name)
            : Strategy(forwarder)
            , ProcessNackTraits(this)
            , m_retxSuppression(RETX_SUPPRESSION_INITIAL, RetxSuppressionExponential::DEFAULT_MULTIPLIER, RETX_SUPPRESSION_MAX)
            , forwarderAcc(forwarder)
        {
            ParsedInstanceName parsed = parseInstanceName(name);
            if (!parsed.parameters.empty()) {
                NDN_THROW(std::invalid_argument("BestRouteStrategy2WithOutFaceWnd does not accept parameters"));
            }
            if (parsed.version && *parsed.version != getStrategyName()[-1].toVersion()) {
                NDN_THROW(std::invalid_argument("BestRouteStrategy2WithOutFaceWnd does not support version "
                                                + to_string(*parsed.version)));
            }
            this->setInstanceName(makeInstanceName(name, getStrategyName()));

            faceCwndMax[0] = 238;    //-/ustc/1&2,258%3=0
            faceCwndMax[1] = 64;
            faceCwndMax[2] = 220;    //-/ustc/0,257%3=2
            for (int i = 0;i < 3;i++) {
                facePending[i] = 0;
                faceDataNum[i] = 0;
                faceFirst[i] = true;
                faceRttEstm[i] = 0;
            }
            forwarder.beforeExpirePendingInterest.connect(
                [this](const nfd::pit::Entry& entry) {
                    for (const pit::OutRecord& i : entry.getOutRecords()) {
                        int index = i.getFace().getId() % 3;
                        this->facePending[index]--;
                        //std::cout << "OnExpire Decrease, face id(" << i.getFace().getId() << "), Pending: " << facePending[i.getFace().getId() % 16 - 2] << std::endl;
                    }
                });
        }

        const Name& BestRouteStrategy2WithECP::getStrategyName()
        {
            static Name strategyName("/localhost/nfd/strategy/best-route2-conges-ecp/%FD%01");
            return strategyName;
        }



        void BestRouteStrategy2WithECP::afterReceiveInterest(const FaceEndpoint& ingress, const Interest& interest,
                                                             const shared_ptr<pit::Entry>& pitEntry)
        {
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
                NFD_LOG_DEBUG(interest << " from=" << ingress << " newPitEntry-to=" << egress);
                int ret;
                if (ret = canSendInt(pitEntry, egress, interest)) {
                    lp::NackHeader nackHeader;
                    if (ret == 1)
                        nackHeader.setReason(lp::NackReason::BUSY);
                    if (ret == 2)
                        nackHeader.setReason(lp::NackReason::FREE);
                    this->sendNack(pitEntry, ingress, nackHeader);
                    this->sendInterest(pitEntry, egress, interest);
                }
                else {
                    lp::NackHeader nackHeader;
                    nackHeader.setReason(lp::NackReason::CONGESTION);
                    this->sendNack(pitEntry, ingress, nackHeader);
                    if (pitEntry->getOutRecords().size() == 0)
                        this->rejectPendingInterest(pitEntry);
                    return;
                }
                return;
            }

            // find an unused upstream with lowest cost except downstream
            it = std::find_if(nexthops.begin(), nexthops.end(), [&](const auto& nexthop) {
                return isNextHopEligible(ingress.face, interest, nexthop, pitEntry, true, time::steady_clock::now());
                              });

            if (it != nexthops.end()) {
                auto egress = FaceEndpoint(it->getFace(), 0);
                int ret;
                if (ret = canSendInt(pitEntry, egress, interest)) {
                    lp::NackHeader nackHeader;
                    if (ret == 1)
                        nackHeader.setReason(lp::NackReason::BUSY);
                    if (ret == 2)
                        nackHeader.setReason(lp::NackReason::FREE);
                    this->sendNack(pitEntry, ingress, nackHeader);
                    this->sendInterest(pitEntry, egress, interest);
                }
                else {
                    lp::NackHeader nackHeader;
                    nackHeader.setReason(lp::NackReason::CONGESTION);
                    this->sendNack(pitEntry, ingress, nackHeader);
                    if (pitEntry->getOutRecords().size() == 0)
                        this->rejectPendingInterest(pitEntry);
                    return;
                }
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
                int ret;
                if (ret = canSendInt(pitEntry, egress, interest)) {
                    lp::NackHeader nackHeader;
                    if (ret == 1)
                        nackHeader.setReason(lp::NackReason::BUSY);
                    if (ret == 2)
                        nackHeader.setReason(lp::NackReason::FREE);
                    this->sendNack(pitEntry, ingress, nackHeader);
                    this->sendInterest(pitEntry, egress, interest);
                }
                else {
                    lp::NackHeader nackHeader;
                    nackHeader.setReason(lp::NackReason::CONGESTION);
                    this->sendNack(pitEntry, ingress, nackHeader);
                    if (pitEntry->getOutRecords().size() == 0)
                        this->rejectPendingInterest(pitEntry);
                    return;
                }
                NFD_LOG_DEBUG(interest << " from=" << ingress << " retransmit-retry-to=" << egress);
            }
        }

        int BestRouteStrategy2WithECP::canSendInt(const shared_ptr<pit::Entry>& pitEntry,
                                                  const FaceEndpoint& egress, const Interest& interest)
        {
            int index = egress.face.getId() % 3;
            int ret = 0;
            if (facePending[index] >= faceCwndMax[index])
                return ret;
            else if (facePending[index] >= 0.5 * faceCwndMax[index])
                ret = 1;
            else
                ret = 2;

            if (pitEntry->getOutRecords().size() == 0)
                facePending[index]++;
            else {
                auto i = pitEntry->getOutRecord(egress.face);
                if (i == pitEntry->out_end())
                    facePending[index]++;
            }
            //std::cout << "OnInterest Increase, " << "face: " << egress << "\tPending: " << facePending[index] << std::endl;
            return ret;
        }

        void BestRouteStrategy2WithECP::afterReceiveNack(const FaceEndpoint& ingress, const lp::Nack& nack,
                                                         const shared_ptr<pit::Entry>& pitEntry)
        {
            int index = ingress.face.getId() % 3;
            if (nack.getReason() != lp::NackReason::FREE && nack.getReason() != lp::NackReason::BUSY) {
                facePending[index]--;
                //std::cout << "OnNack Decrease, " << "face: " << ingress << "\tPending: " << facePending[index] << std::endl;
            }
            this->processNack(ingress.face, nack, pitEntry);
        }

        void BestRouteStrategy2WithECP::beforeSatisfyInterest(const shared_ptr<pit::Entry>& pitEntry,
                                                              const FaceEndpoint& ingress, const Data& data)
        {
            if (ingress.face.getId() != face::FACEID_CONTENT_STORE) {
                int index = ingress.face.getId() % 3;
                facePending[index]--;
                faceDataNum[index] = ingress.face.getCounters().nInData;
                //std::cout << "OnData Decrease, " << "face: " << ingress
                //    << "\tPending: " << facePending[index]
                //    << "\tMeanRtt: " << faceRttEstm[index]
                //    << "\tDataNum: " << faceDataNum[index]
                //    << std::endl;
            }
            Strategy::beforeSatisfyInterest(pitEntry, ingress, data);
        }

        void BestRouteStrategy2WithECP::afterContentStoreHit(const shared_ptr<pit::Entry>& pitEntry,
                                                             const FaceEndpoint& ingress, const Data& data)
        {
            //int index = ingress.face.getId() % 16 - 2;
            //facePending[index]++;
            Strategy::afterContentStoreHit(pitEntry, ingress, data);
        }

    } // namespace fw
} // namespace nfd
