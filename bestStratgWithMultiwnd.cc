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

#include "bestStratgWithMultiwnd.h"
#include "ns3/ndnSIM/NFD/daemon/fw/algorithm.hpp"
#include "ns3/ndnSIM/NFD/daemon/common/logger.hpp"
#include "ns3/ndnSIM/NFD/core/common.hpp"

namespace nfd {
    namespace fw {

        NFD_LOG_INIT(BestRouteStrategy2WithOutFaceWnd);
        NFD_REGISTER_STRATEGY(BestRouteStrategy2WithOutFaceWnd);

        const time::milliseconds BestRouteStrategy2WithOutFaceWnd::RETX_SUPPRESSION_INITIAL(10);
        const time::milliseconds BestRouteStrategy2WithOutFaceWnd::RETX_SUPPRESSION_MAX(250);

        BestRouteStrategy2WithOutFaceWnd::BestRouteStrategy2WithOutFaceWnd(Forwarder& forwarder, const Name& name)
            : Strategy(forwarder)
            , ProcessNackTraits(this)
            , m_retxSuppression(RETX_SUPPRESSION_INITIAL, RetxSuppressionExponential::DEFAULT_MULTIPLIER, RETX_SUPPRESSION_MAX)
            , forwarderAcc(forwarder), dTimer(ns3::Seconds(0.1), 1024, this)
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

            for (int i = 0;i < 3;i++) {
                faceCwndMax[i] = 64;
                facePending[i] = 0;
                faceDataNum[i] = 0;
                faceFirst[i] = true;
                faceRttEstm[i] = 0;
            }
            forwarder.beforeExpirePendingInterest.connect(
                [this](const nfd::pit::Entry& entry) {
                    for (const pit::OutRecord& i : entry.getOutRecords()) {
                        int index = i.getFace().getId() % 16 - 2;
                        this->facePending[index]--;
                        this->faceRttRec[index].erase(entry.getName());
                        //std::cout << "OnExpire Decrease, face id(" << i.getFace().getId() << "), Pending: " << facePending[i.getFace().getId() % 16 - 2] << std::endl;
                    }
                });
        }

        const Name& BestRouteStrategy2WithOutFaceWnd::getStrategyName()
        {
            static Name strategyName("/localhost/nfd/strategy/best-route2-conges-multiport/%FD%01");
            return strategyName;
        }



        void BestRouteStrategy2WithOutFaceWnd::afterReceiveInterest(const FaceEndpoint& ingress, const Interest& interest,
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
                if (canSendInt(pitEntry, egress, interest))
                    this->sendInterest(pitEntry, egress, interest);
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
                if (canSendInt(pitEntry, egress, interest))
                    this->sendInterest(pitEntry, egress, interest);
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
                if (canSendInt(pitEntry, egress, interest))
                    this->sendInterest(pitEntry, egress, interest);
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

        bool BestRouteStrategy2WithOutFaceWnd::canSendInt(const shared_ptr<pit::Entry>& pitEntry,
                                                          const FaceEndpoint& egress, const Interest& interest)
        {
            int index = egress.face.getId() % 16 - 2;
            if (facePending[index] >= faceCwndMax[index])
                return false;

            if (pitEntry->getOutRecords().size() == 0)
                facePending[index]++;
            else {
                auto i = pitEntry->getOutRecord(egress.face);
                if (i == pitEntry->out_end())
                    facePending[index]++;
            }
            faceRttRec[index].insert(std::make_pair(interest.getName(), ns3::Simulator::Now()));

            //std::cout << "OnInterest Increase, " << "face: " << egress << "\tPending: " << facePending[index] << std::endl;
            return true;
        }

        void BestRouteStrategy2WithOutFaceWnd::afterReceiveNack(const FaceEndpoint& ingress, const lp::Nack& nack,
                                                                const shared_ptr<pit::Entry>& pitEntry)
        {
            int index = ingress.face.getId() % 16 - 2;
            if (nack.getReason() != lp::NackReason::FREE && nack.getReason() != lp::NackReason::BUSY) {
                facePending[index]--;
                faceRttRec[index].erase(nack.getInterest().getName());
                //std::cout << "OnNack Decrease, " << "face: " << ingress << "\tPending: " << facePending[index] << std::endl;
            }
            this->processNack(ingress.face, nack, pitEntry);
        }

        void BestRouteStrategy2WithOutFaceWnd::beforeSatisfyInterest(const shared_ptr<pit::Entry>& pitEntry,
                                                                     const FaceEndpoint& ingress, const Data& data)
        {
            if (ingress.face.getId() != face::FACEID_CONTENT_STORE) {
                int index = ingress.face.getId() % 16 - 2;
                facePending[index]--;
                faceDataNum[index] = ingress.face.getCounters().nInData;
                ns3::Time difftime = ns3::Simulator::Now() - faceRttRec[index][data.getName()];
                if (faceFirst[index]) {
                    faceRttEstm[index] = difftime.GetSeconds();
                    faceFirst[index] = false;
                }
                if (difftime.GetSeconds() != 0) {
                    faceRttEstm[index] = faceRttEstm[index] * 0.75 + 0.25 * difftime.GetSeconds();
                    faceRttRec[index].erase(data.getName());
                }

                //std::cout << "OnData Decrease, " << "face: " << ingress
                //    << "\tPending: " << facePending[index]
                //    << "\tMeanRtt: " << faceRttEstm[index]
                //    << "\tDataNum: " << faceDataNum[index]
                //    << std::endl;
            }
            Strategy::beforeSatisfyInterest(pitEntry, ingress, data);
        }

        void BestRouteStrategy2WithOutFaceWnd::afterContentStoreHit(const shared_ptr<pit::Entry>& pitEntry,
                                                                    const FaceEndpoint& ingress, const Data& data)
        {
            //int index = ingress.face.getId() % 16 - 2;
            //facePending[index]++;
            Strategy::afterContentStoreHit(pitEntry, ingress, data);
        }

        ddpgTimerTrans::ddpgTimerTrans(ns3::Time interval, uint16_t id,
                                       BestRouteStrategy2WithOutFaceWnd* ptr)
            :ns3::Ns3AIRL<edgeDDPGEnv, edgeDDPGAct>(id)
            , intervalt(interval)
        {
            Ping(intervalt);
            SetFunction(BestRouteStrategy2WithOutFaceWnd::callback);
            //SetArguments<BestRouteStrategy2WithOutFaceWnd*>(ptr);
            SetCond(2, 0);
        }

        void BestRouteStrategy2WithOutFaceWnd::callback(BestRouteStrategy2WithOutFaceWnd* ptr)
        {
            std::cout << std::setw(64) << std::setfill('=') << '=' << std::setfill(' ') << std::endl;
            std::cout << "time:" << ns3::Simulator::Now().GetSeconds() << std::endl;
            std::cout << std::setiosflags(std::ios::left)
                << std::setw(12) << "face"
                << std::setw(12) << "Cwndmax"
                << std::setw(12) << "Pending"
                << std::setw(12) << "DataNum"
                << std::setw(12) << "RttEstm"
                << std::endl;
            for (int i = 0;i < 3;i++) {
                std::cout << std::setiosflags(std::ios::left)
                    << std::setw(12) << i
                    << std::setw(12) << ptr->faceCwndMax[i]
                    << std::setw(12) << ptr->facePending[i]
                    << std::setw(12) << ptr->faceDataNum[i]
                    << std::setw(12) << ptr->faceRttEstm[i]
                    << std::endl;
            }
            std::cout << std::setw(64) << std::setfill('-') << '-' << std::setfill(' ') << std::endl;
            ptr->dTimer.env.time_s = ns3::Simulator::Now().GetSeconds();
            memcpy(ptr->dTimer.env.faceCwndMax, ptr->faceCwndMax, sizeof(faceCwndMax));
            memcpy(ptr->dTimer.env.facePending, ptr->facePending, sizeof(facePending));
            memcpy(ptr->dTimer.env.faceDataNum, ptr->faceDataNum, sizeof(faceDataNum));
            memcpy(ptr->dTimer.env.faceRttEstm, ptr->faceRttEstm, sizeof(faceRttEstm));
            ptr->dTimer.sendEnv();
            ptr->dTimer.getAct();
            for (int i : {0, 1, 2}) {
                ptr->faceCwndMax[i] *= ptr->dTimer.act.faceActVal[i];
                if (ptr->faceCwndMax[i] < 1)
                    ptr->faceCwndMax[i] = 1;
                else if (ptr->faceCwndMax[i] > 500)
                    ptr->faceCwndMax[i] = 500;
            }
            std::cout << std::setiosflags(std::ios::left)
                << std::setw(12) << "face"
                << std::setw(12) << "Val"
                << std::endl;
            for (int i = 0;i < 3;i++) {
                std::cout << std::setiosflags(std::ios::left)
                    << std::setw(12) << i
                    << std::setw(12) << ptr->faceCwndMax[i]
                    << std::endl;
            }
            std::cout << std::setw(64) << std::setfill('=') << '=' << std::setfill(' ') << std::endl;
            ptr->dTimer.Ping(ptr->dTimer.intervalt);
        }

        void ddpgTimerTrans::sendEnv()
        {
            auto shmEnv = EnvSetterCond();
            memcpy(shmEnv, &env, sizeof(edgeDDPGEnv));
            SetCompleted();
        }

        void ddpgTimerTrans::getAct()
        {
            auto shmAct = ActionGetterCond();
            memcpy(&act, shmAct, sizeof(edgeDDPGAct));
            GetCompleted();
        }
    } // namespace fw
} // namespace nfd
