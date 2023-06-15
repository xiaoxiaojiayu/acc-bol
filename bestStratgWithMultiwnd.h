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

#ifndef NFD_DAEMON_FW_BEST_ROUTE_STRATEGY2_WITH_CONGES_HPP
#define NFD_DAEMON_FW_BEST_ROUTE_STRATEGY2_WITH_CONGES_HPP

#include "ns3/ndnSIM/NFD/daemon/fw/strategy.hpp"
#include "ns3/ndnSIM/NFD/daemon/fw/process-nack-traits.hpp"
#include "ns3/ndnSIM/NFD/daemon/fw/retx-suppression-exponential.hpp"
#include "ns3/watchdog.h"
#include "ns3/ns3-ai-module.h"

namespace nfd {
    namespace fw {

        class BestRouteStrategy2WithOutFaceWnd;

        typedef struct {
            double time_s;
            double faceCwndMax[3];
            int facePending[3];
            int faceDataNum[3];
            double faceRttEstm[3];
        }Packed edgeDDPGEnv;

        typedef struct {
            double faceActVal[3];
        }Packed edgeDDPGAct;

        class ddpgTimerTrans
            : ns3::Watchdog
            , ns3::Ns3AIRL<edgeDDPGEnv, edgeDDPGAct> {
        public:
            ddpgTimerTrans(ns3::Time interval, uint16_t id, BestRouteStrategy2WithOutFaceWnd* ptr);
            //~ddpgTimerTrans();
            void sendEnv();
            void getAct();
            edgeDDPGEnv env;
            edgeDDPGAct act;
        private:
            const ns3::Time intervalt;
            friend BestRouteStrategy2WithOutFaceWnd;
        };


        class BestRouteStrategy2WithOutFaceWnd
            : public Strategy
            , public ProcessNackTraits<BestRouteStrategy2WithOutFaceWnd> {
        public:
            explicit
                BestRouteStrategy2WithOutFaceWnd(Forwarder& forwarder, const Name& name = getStrategyName());

            static const Name&
                getStrategyName();

            bool canSendInt(const shared_ptr<pit::Entry>& pitEntry, const FaceEndpoint& egress, const Interest& interest);

            void afterReceiveInterest(const FaceEndpoint& ingress, const Interest& interest,
                                      const shared_ptr<pit::Entry>& pitEntry) override;

            void afterReceiveNack(const FaceEndpoint& ingress, const lp::Nack& nack,
                                  const shared_ptr<pit::Entry>& pitEntry)override;

            void beforeSatisfyInterest(const shared_ptr<pit::Entry>& pitEntry,
                                       const FaceEndpoint& ingress, const Data& data);

            void afterContentStoreHit(const shared_ptr<pit::Entry>& pitEntry,
                                      const FaceEndpoint& ingress, const Data& data);

            static void callback(BestRouteStrategy2WithOutFaceWnd* ptr);

        private:
            Forwarder& forwarderAcc;
            ddpgTimerTrans dTimer;
            double faceCwndMax[3];
            int facePending[3];
            int faceDataNum[3];
            double faceRttEstm[3];
            bool faceFirst[3];
            std::map<Name, ns3::Time> faceRttRec[3];
            

            static const time::milliseconds RETX_SUPPRESSION_INITIAL;
            static const time::milliseconds RETX_SUPPRESSION_MAX;
            RetxSuppressionExponential m_retxSuppression;

            friend ProcessNackTraits<BestRouteStrategy2WithOutFaceWnd>;
        };

    } 
} 
#endif 
