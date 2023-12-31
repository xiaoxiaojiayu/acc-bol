/* -*- Mode:C++; c-file-style:"gnu"; indent-tabs-mode:nil; -*- */
/*
 * Copyright (c) 2011-2018  Regents of the University of California.
 *
 * This file is part of ndnSIM. See AUTHORS for complete list of ndnSIM authors and
 * contributors.
 *
 * ndnSIM is free software: you can redistribute it and/or modify it under the terms
 * of the GNU General Public License as published by the Free Software Foundation,
 * either version 3 of the License, or (at your option) any later version.
 *
 * ndnSIM is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 * PURPOSE.  See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * ndnSIM, e.g., in COPYING.md file.  If not, see <http://www.gnu.org/licenses/>.
 **/

#include "pcon1.hpp"
#include <iomanip>
#include <iostream>
#include <limits>
#include <queue>
#include <random>
#include <unordered_set>

#include "gsl/gsl_fit.h"
#include "ns3/assert.h"
#include "ns3/callback.h"
#include "ns3/core-module.h"
#include "ns3/double.h"
#include "ns3/log.h"
#include "ns3/scheduler.h"
#include "ns3/string.h"
#include "ns3/uinteger.h"

NS_LOG_COMPONENT_DEFINE("ndn.ConsumerPcon1");

namespace ns3 {
namespace ndn {

NS_OBJECT_ENSURE_REGISTERED(ConsumerPcon1);

constexpr double ConsumerPcon1::CUBIC_C;
constexpr uint32_t ConsumerPcon1::BIC_MAX_INCREMENT;
constexpr uint32_t ConsumerPcon1::BIC_LOW_WINDOW;

TypeId
ConsumerPcon1::GetTypeId()
{
  static TypeId tid =
    TypeId("ns3::ndn::ConsumerPcon1")
      .SetGroupName("Ndn")
      .SetParent<ConsumerWindow>()
      .AddConstructor<ConsumerPcon1>()

      .AddAttribute("CcAlgorithm",
                    "Specify which window adaptation algorithm to use (AIMD, BIC, or CUBIC)",
                    EnumValue(CcAlgorithm::AIMD),
                    MakeEnumAccessor(&ConsumerPcon1::m_ccAlgorithm),
                    MakeEnumChecker(CcAlgorithm::AIMD, "AIMD", CcAlgorithm::BIC, "BIC",
                                    CcAlgorithm::CUBIC, "CUBIC"))

      .AddAttribute("Beta",
                    "TCP Multiplicative Decrease factor",
                    DoubleValue(0.5),
                    MakeDoubleAccessor(&ConsumerPcon1::m_beta),
                    MakeDoubleChecker<double>())

      .AddAttribute("CubicBeta",
                    "TCP CUBIC Multiplicative Decrease factor",
                    DoubleValue(0.8),
                    MakeDoubleAccessor(&ConsumerPcon1::m_cubicBeta),
                    MakeDoubleChecker<double>())

      .AddAttribute("AddRttSuppress",
                    "Minimum number of RTTs (1 + this factor) between window decreases",
                    DoubleValue(0.5), // This default value was chosen after manual testing
                    MakeDoubleAccessor(&ConsumerPcon1::m_addRttSuppress),
                    MakeDoubleChecker<double>())

      .AddAttribute("ReactToCongestionMarks",
                    "If true, process received congestion marks",
                    BooleanValue(true),
                    MakeBooleanAccessor(&ConsumerPcon1::m_reactToCongestionMarks),
                    MakeBooleanChecker())

      .AddAttribute("UseCwa",
                    "If true, use Conservative Window Adaptation",
                    BooleanValue(true),
                    MakeBooleanAccessor(&ConsumerPcon1::m_useCwa),
                    MakeBooleanChecker())

      .AddAttribute("UseCubicFastConvergence",
                    "If true, use TCP CUBIC Fast Convergence",
                    BooleanValue(false),
                    MakeBooleanAccessor(&ConsumerPcon1::m_useCubicFastConv),
                    MakeBooleanChecker())
                    
      .AddAttribute("RandomPrefix", "",
                StringValue(""),
                MakeStringAccessor(&ConsumerPcon1::SetRandomPrefix),
                MakeStringChecker());    ;

  return tid;
}

ConsumerPcon1::ConsumerPcon1()
  : m_ssthresh(8000)
  , m_highData(0)
  , m_recPoint(0.0)
  , m_cubicWmax(0)
  , m_cubicLastWmax(0)
  , m_cubicLastDecrease(time::steady_clock::now())
  , m_bicMinWin(0)
  , m_bicMaxWin(std::numeric_limits<double>::max())
  , m_bicTargetWin(0)
  , m_bicSsCwnd(0)
  , m_bicSsTarget(0)
  , m_flag(0)
  , m_isBicSs(false)
{
}

void
ConsumerPcon1::OnData(shared_ptr<const Data> data)
{
  Consumer::OnData(data);

  uint64_t sequenceNum = data->getName().get(-1).toSequenceNumber();

  // Set highest received Data to sequence number
  if (m_highData < sequenceNum) {
    m_highData = sequenceNum;
  }

  if (data->getCongestionMark() > 23) {
    if (m_reactToCongestionMarks) {
      NS_LOG_DEBUG("Received congestion mark: " << data->getCongestionMark());
      WindowDecrease();
    }
    else {
      NS_LOG_DEBUG("Ignored received congestion mark: " << data->getCongestionMark());
    }
  }
  else {
    WindowIncrease();
  }

  if (m_inFlight > static_cast<uint32_t>(0)) {
    m_inFlight--;
  }

  NS_LOG_DEBUG("Window: " << m_window << ", InFlight: " << m_inFlight);

  ScheduleNextPacket();
}

void
ConsumerPcon1::OnTimeout(uint32_t sequenceNum)
{
  m_flag =1;
  WindowDecrease();

  if (m_inFlight > static_cast<uint32_t>(0)) {
    m_inFlight--;
  }

  NS_LOG_DEBUG("Window: " << m_window << ", InFlight: " << m_inFlight);

  Consumer::OnTimeout(sequenceNum);
}

void
ConsumerPcon1::WindowIncrease()
{
  if (m_ccAlgorithm == CcAlgorithm::AIMD) {
    if (m_window < m_ssthresh) {
      m_window += 1.0;
    }
    else {
      m_window += (1.0 / m_window);
    }
  }
  else if (m_ccAlgorithm == CcAlgorithm::CUBIC) {
    CubicIncrease();
  }
  else if (m_ccAlgorithm == CcAlgorithm::BIC) {
    BicIncrease();
  }
  else {
    BOOST_ASSERT_MSG(false, "Unknown CC Algorithm");
  }
  NS_LOG_DEBUG("Window size increased to " << m_window);
}

void
ConsumerPcon1::WindowDecrease()
{
  if (!m_useCwa || m_highData > m_recPoint) {
    const double diff = m_seq - m_highData;
    BOOST_ASSERT(diff > 0);

    m_recPoint = m_seq + (m_addRttSuppress * diff);

    if (m_ccAlgorithm == CcAlgorithm::AIMD) {
      // Normal TCP Decrease:
      m_ssthresh = m_window * m_beta;
      m_window = m_ssthresh;
    }
    else if (m_ccAlgorithm == CcAlgorithm::CUBIC) {
      CubicDecrease();
    }
    else if (m_ccAlgorithm == CcAlgorithm::BIC) {
      BicDecrease();
    }
    else {
      BOOST_ASSERT_MSG(false, "Unknown CC Algorithm");
    }

    // Window size cannot be reduced below initial size
    if (m_window < m_initialWindow) {
      m_window = m_initialWindow;
    }

    NS_LOG_DEBUG("Window size decreased to " << m_window);
  }
  else {
    NS_LOG_DEBUG("Window decrease suppressed, HighData: " << m_highData << ", RecPoint: " << m_recPoint);
  }
}


void
ConsumerPcon1::BicIncrease()
{
  if (m_window < BIC_LOW_WINDOW) {
    // Normal TCP AIMD behavior
    if (m_window < m_ssthresh) {
      m_window = m_window + 1;
    }
    else {
      m_window = m_window + 1.0 / m_window;
    }
  }
  else if (!m_isBicSs) {
    // Binary increase
    if (m_bicTargetWin - m_window < BIC_MAX_INCREMENT) { // Binary search
      m_window += (m_bicTargetWin - m_window) / m_window;
    }
    else {
      m_window += BIC_MAX_INCREMENT / m_window; // Additive increase
    }
    // FIX for equal double values.
    if (m_window + 0.00001 < m_bicMaxWin) {
      m_bicMinWin = m_window;
      m_bicTargetWin = (m_bicMaxWin + m_bicMinWin) / 2;
    }
    else {
      m_isBicSs = true;
      m_bicSsCwnd = 1;
      m_bicSsTarget = m_window + 1.0;
      m_bicMaxWin = std::numeric_limits<double>::max();
    }
  }
  else {
    // BIC slow start
    m_window += m_bicSsCwnd / m_window;
    if (m_window >= m_bicSsTarget) {
      m_bicSsCwnd = 2 * m_bicSsCwnd;
      m_bicSsTarget = m_window + m_bicSsCwnd;
    }
    if (m_bicSsCwnd >= BIC_MAX_INCREMENT) {
      m_isBicSs = false;
    }
  }
}

void
ConsumerPcon1::BicDecrease()
{
  // BIC Decrease
  if (m_window >= BIC_LOW_WINDOW) {
    auto prev_max = m_bicMaxWin;
    m_bicMaxWin = m_window;
    m_window = m_window * m_cubicBeta;
    m_bicMinWin = m_window;
    if (prev_max > m_bicMaxWin) {
      // Fast Convergence
      m_bicMaxWin = (m_bicMaxWin + m_bicMinWin) / 2;
    }
    m_bicTargetWin = (m_bicMaxWin + m_bicMinWin) / 2;
  }
  else {
    // Normal TCP Decrease:
    m_ssthresh = m_window * m_cubicBeta;
    m_window = m_ssthresh;
  }
}


void
ConsumerPcon1::CubicIncrease()
{
  // 1. Time since last congestion event in Seconds
  const double t = time::duration_cast<time::microseconds>(
                     time::steady_clock::now() - m_cubicLastDecrease).count() / 1e6;

  // 2. Time it takes to increase the window to cubic_wmax
  // K = cubic_root(W_max*(1-beta_cubic)/C) (Eq. 2)
  const double k = std::cbrt(m_cubicWmax * (1 - m_cubicBeta) / CUBIC_C);

  // 3. Target: W_cubic(t) = C*(t-K)^3 + W_max (Eq. 1)
  const double w_cubic = CUBIC_C * std::pow(t - k, 3) + m_cubicWmax;

  // 4. Estimate of Reno Increase (Currently Disabled)
  //  const double rtt = m_rtt->GetCurrentEstimate().GetSeconds();
  //  const double w_est = m_cubic_wmax*m_beta + (3*(1-m_beta)/(1+m_beta)) * (t/rtt);
  constexpr double w_est = 0.0;

  // Actual adaptation
  if (m_window < m_ssthresh) {
    m_window += 1.0;
  }
  else {
    BOOST_ASSERT(m_cubicWmax > 0);

    double cubic_increment = std::max(w_cubic, w_est) - m_window;
    // Cubic increment must be positive:
    // Note: This change is not part of the RFC, but I added it to improve performance.
    if (cubic_increment < 0) {
      cubic_increment = 0.0;
    }
    m_window += cubic_increment / m_window;
  }
}


void
ConsumerPcon1::CubicDecrease()
{
  // This implementation is ported from https://datatracker.ietf.org/doc/rfc8312/

  const double FAST_CONV_DIFF = 1.0; // In percent

  // A flow remembers the last value of W_max,
  // before it updates W_max for the current congestion event.

  // Current w_max < last_wmax
  if (m_useCubicFastConv && m_window < m_cubicLastWmax * (1 - FAST_CONV_DIFF / 100)) {
    m_cubicLastWmax = m_window;
    m_cubicWmax = m_window * (1.0 + m_cubicBeta) / 2.0;
  }
  else {
    // Save old cwnd as w_max:
    m_cubicLastWmax = m_window;
    m_cubicWmax = m_window;
  }
  if(m_flag==0)
  {
      m_ssthresh = m_window * m_cubicBeta;
       m_ssthresh = std::max<double>(m_ssthresh, m_initialWindow);
      m_window = m_ssthresh;

      m_cubicLastDecrease = time::steady_clock::now();
  }
  else{
    m_ssthresh = m_window * 0.5;
       m_ssthresh = std::max<double>(m_ssthresh, m_initialWindow);
      m_window = m_ssthresh;

      m_cubicLastDecrease = time::steady_clock::now();
      m_flag = 0;
  }
  
}


void ConsumerPcon1::SetRandomPrefix(std::string random_prefix)  {
       if (!random_prefix.empty())
       {
           double acml = 0;
           std::size_t bpos = 0, epos = 0;
           do
           {
                epos = random_prefix.find(',', bpos);
                std::string kvpair = random_prefix.substr(bpos, epos - bpos);
                bpos = epos + 1;
                acml += std::stod(kvpair.substr(kvpair.find('=') + 1));
                intervalp.push_back(acml);
                intervalmap.emplace(std::make_pair(acml, kvpair.substr(0, kvpair.find('='))));
            } while (bpos != std::string::npos + 1);
             randprefix = true;
       }
        else
            randprefix = false;
}



void ConsumerPcon1::ScheduleNextPacket()  {
      if (randprefix)
       {
           static std::random_device rdev;
           static std::mt19937 reng(rdev());
           static std::uniform_real_distribution<> u(0, intervalp.back());
           double rnum = u(reng);
           int i;
           for (i = 0; i < intervalp.size(); i++)
                if (rnum <= intervalp[i])
                   break;
          this->SetAttribute("Prefix", StringValue(intervalmap.find(intervalp[i])->second));
       }
	
	ConsumerWindow::ScheduleNextPacket();
}

// void ConsumerPcon1::changeprefix() {
//     NS_LOG_DEBUG(m_interestName);
//     if (m_interestName == "/ustc")
//         m_interestName = "/ustc1";
//     else
//         m_interestName = "/ustc";
//     NS_LOG_DEBUG(m_interestName);
// }

} // namespace ndn
} // namespace ns3
