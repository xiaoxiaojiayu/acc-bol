/* -*- Mode:C++; c-file-style:"gnu"; indent-tabs-mode:nil; -*- */
/*
 * Copyright (c) 2011-2018  Regents of the University of California.
 *
 * This file is part of ndnSIM. See AUTHORS for complete list of ndnSIM authors
 *and contributors.
 *
 * ndnSIM is free software: you can redistribute it and/or modify it under the
 *terms of the GNU General Public License as published by the Free Software
 *Foundation, either version 3 of the License, or (at your option) any later
 *version.
 *
 * ndnSIM is distributed in the hope that it will be useful, but WITHOUT ANY
 *WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
 *A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * ndnSIM, e.g., in COPYING.md file.  If not, see
 *<http://www.gnu.org/licenses/>.
 **/

#include "ndn-wifi-net-device-transport.hpp"

#include <ns3/ndnSIM/ndn-cxx/data.hpp>
#include <ns3/ndnSIM/ndn-cxx/encoding/block.hpp>
#include <ns3/ndnSIM/ndn-cxx/interest.hpp>

#include "ns3/ap-wifi-mac.h"
#include "ns3/ndnSIM//utils/ndn-ns3-packet-tag.hpp"
#include "ns3/ndnSIM/helper/ndn-stack-helper.hpp"
#include "ns3/ndnSIM/model/ndn-block-header.hpp"
#include "ns3/queue.h"
#include "ns3/sta-wifi-mac.h"
#include "ns3/wifi-mac-queue.h"
#include "ns3/wifi-mac.h"
#include "ns3/yans-wifi-channel.h"

NS_LOG_COMPONENT_DEFINE("ndn.WifiNetDeviceTransport");

namespace ns3 {
namespace ndn {

WifiNetDeviceTransport::WifiNetDeviceTransport(
    Ptr<Node> node, const Ptr<WifiNetDevice>& netDevice,
    const std::string& localUri, const std::string& remoteUri,
    ::ndn::nfd::FaceScope scope, ::ndn::nfd::FacePersistency persistency,
    ::ndn::nfd::LinkType linkType)
    : m_netDevice(netDevice), m_node(node), macqueue(nullptr) {
    this->setLocalUri(FaceUri(localUri));
    this->setRemoteUri(FaceUri(remoteUri));
    this->setScope(scope);
    this->setPersistency(persistency);
    this->setLinkType(linkType);
    this->setMtu(m_netDevice->GetMtu());  // Use the MTU of the netDevice

    local_addr = getLocalUri().getHost().c_str();
    remote_addr = getRemoteUri().getHost().c_str();

    // Get send queue capacity for congestion marking
    PointerValue txQueueAttribute;
    if (m_netDevice->GetAttributeFailSafe("TxQueue", txQueueAttribute)) {
        NS_LOG_DEBUG("get wifi net device queue");
        Ptr<ns3::QueueBase> txQueue = txQueueAttribute.Get<ns3::QueueBase>();
        // must be put into bytes mode queue

        auto size = txQueue->GetMaxSize();
        NS_LOG_DEBUG("get wifi netdevice queue, size: " << size);
        if (size.GetUnit() == BYTES) {
            this->setSendQueueCapacity(size.GetValue());
        } else {
            // don't know the exact size in bytes, guessing based on "standard"
            // packet size
            this->setSendQueueCapacity(size.GetValue() * 1500);
        }
    }

    auto mac = m_netDevice->GetMac();
    if (mac->GetTypeOfStation() == ns3::TypeOfStation::AP) {
        macqueue = DynamicCast<ns3::ApWifiMac>(mac)->GetTxopQueue(
            ns3::AcIndex::AC_BE);  // ns3::AcIndex::AC_BE// (ns3::AcIndex)i
    }

    NS_LOG_FUNCTION(
        this << "Creating an ndnSIM transport instance for netDevice with URI"
             << this->getLocalUri());

    NS_ASSERT_MSG(m_netDevice != 0,
                  "NetDeviceFace needs to be assigned a valid NetDevice");

    m_node->RegisterProtocolHandler(
        MakeCallback(&WifiNetDeviceTransport::receiveFromNetDevice, this),
        L3Protocol::ETHERNET_FRAME_TYPE, m_netDevice,
        true /*promiscuous mode*/);
}

WifiNetDeviceTransport::~WifiNetDeviceTransport() { NS_LOG_FUNCTION_NOARGS(); }

ssize_t WifiNetDeviceTransport::getSendQueueLength() {
    // PointerValue txQueueAttribute;
    // if (m_netDevice->GetAttributeFailSafe("TxQueue", txQueueAttribute)) {
    //     Ptr<ns3::QueueBase> txQueue = txQueueAttribute.Get<ns3::QueueBase>();
    //     return txQueue->GetNBytes();
    // } else {
    //     return nfd::face::QUEUE_UNSUPPORTED;
    // }
    ns3::Ptr<ns3::WifiMacQueue> macqueue;
    auto mac = m_netDevice->GetMac();
    if (mac->GetTypeOfStation() == ns3::TypeOfStation::AP) {
        macqueue = DynamicCast<ns3::ApWifiMac>(mac)->GetTxopQueue(
            ns3::AcIndex::AC_BE);  // ns3::AcIndex::AC_BE// (ns3::AcIndex)i
    }
    if (macqueue != nullptr) {
        //std::cout <<macqueue<<" "<< macqueue->GetNPackets() << "\n";
        return macqueue->GetNPackets();
    } else {
        return -1;
    }
    // else if (mac->GetTypeOfStation() == ns3::TypeOfStation::STA) {
    //     // macqueue = DynamicCast<ns3::StaWifiMac>(mac)->GetTxopQueue(
    //     //     ns3::AcIndex::AC_BE);
    //     for (int i = 0; i < 1; i++) {
    //         macqueue[i] = DynamicCast<ns3::StaWifiMac>(mac)->GetTxopQueue(
    //             ns3::AcIndex::AC_BE);  // ns3::AcIndex::AC_BE
    //     }
    // }
    // for (int i = 0; i < 1; i++) {
    //     if (macqueue[i] != nullptr)
    //         // macqueue[i]->GetNPackets();
    //         std::cout << "maxqueue "<<macqueue[i]<<" i: " << i << "
    //         macqueue[i]->GetNPackets(): "
    //              << macqueue[i]->GetNPackets()<<"\n";
    // }
}

void WifiNetDeviceTransport::doClose() {
    NS_LOG_FUNCTION(this << "Closing transport for netDevice with URI"
                         << this->getLocalUri());

    // set the state of the transport to "CLOSED"
    this->setState(nfd::face::TransportState::CLOSED);
}

void WifiNetDeviceTransport::doSend(const Block& packet,
                                    const nfd::EndpointId& endpoint) {
    NS_LOG_FUNCTION(this << "Sending packet from netDevice with URI"
                         << this->getLocalUri());

    // convert NFD packet to NS3 packet
    BlockHeader header(packet);

    Ptr<ns3::Packet> ns3Packet = Create<ns3::Packet>();
    ns3Packet->AddHeader(header);

    // auto ch = DynamicCast<YansWifiChannel>(m_netDevice->GetChannel());

    // m_netDevice->GetMac();

    // send the NS3 packet
    m_netDevice->Send(ns3Packet, remote_addr,  // m_netDevice->GetBroadcast(),
                      L3Protocol::ETHERNET_FRAME_TYPE);
}

// callback
void WifiNetDeviceTransport::receiveFromNetDevice(
    Ptr<NetDevice> device, Ptr<const ns3::Packet> p, uint16_t protocol,
    const Address& from, const Address& to, NetDevice::PacketType packetType) {
    NS_LOG_FUNCTION(device << p << protocol << from << to << packetType);

    // Convert NS3 packet to NFD packet
    Ptr<ns3::Packet> packet = p->Copy();

    BlockHeader header;
    packet->RemoveHeader(header);

    auto mac_from = Mac48Address::ConvertFrom(from);
    bool shouldup = (mac_from == remote_addr);
    NS_LOG_DEBUG(mac_from << " == " << remote_addr << "? it is " << shouldup);

    if (shouldup) {
        this->receive(std::move(header.getBlock()));
    }
}

Ptr<NetDevice> WifiNetDeviceTransport::GetNetDevice() const {
    return m_netDevice;
}

}  // namespace ndn
}  // namespace ns3
