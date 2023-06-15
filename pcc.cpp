#include "pcc.hpp"

#include <queue>
#include <unordered_set>
#include "ns3/callback.h"
#include "gsl/gsl_fit.h"
#include "ns3/assert.h"
#include "ns3/core-module.h"
#include "ns3/log.h"
#include "ns3/scheduler.h"
#include <limits>
#include <iostream>
#include <iomanip>
#include <random>
#include "ns3/string.h"
#include "ns3/uinteger.h"
#include "ns3/double.h"


NS_LOG_COMPONENT_DEFINE("ndn.PCC");

namespace ns3 {
namespace ndn {

NS_OBJECT_ENSURE_REGISTERED(ConsumerPCC);

TypeId ConsumerPCC::GetTypeId() {
    static TypeId tid =
        TypeId("ns3::ndn::ConsumerPCC")
            .SetGroupName("Ndn")
            .SetParent<ConsumerCbr>()
            .AddConstructor<ConsumerPCC>()
            .AddAttribute("rate", "", DoubleValue(100.0),
                          MakeDoubleAccessor(&ConsumerPCC::setrate),
                          MakeDoubleChecker<double>())
            .AddAttribute("epsilon", "", DoubleValue(0.05),
                          MakeDoubleAccessor(&ConsumerPCC::epsilon),
                          MakeDoubleChecker<double>())
	    //.AddAttribute("Frequency", "",
            //                      DoubleValue(100.0),
            //                      MakeDoubleAccessor(&ConsumerPCC::SetFrequency),
            //                      MakeDoubleChecker<double>())
	    //.AddAttribute("Randomize", "",
            //                      StringValue("uniform"),
            //                      MakeStringAccessor(&ConsumerPCC::SetRandomize),
            //                      MakeStringChecker())
	    .AddAttribute("RandomPrefix", "",
                                  StringValue(""),
                                  MakeStringAccessor(&ConsumerPCC::SetRandomPrefix),
                                  MakeStringChecker());         

    return tid;
}

ConsumerPCC::ConsumerPCC()
    : ConsumerCbr(),
      rate(10.0),
      rtt(MilliSeconds(60)),
      MI(MilliSeconds(60)),
      epsilon(0.05),
      theta0(1),
      omega0(0.05),
      delta(0.1),
      omegaboundk(0),
      probestate(0),
      monitorstate(0),
      lastdiffpositive(false),
      contk(1){
    m_firstInterestDataDelay.ConnectWithoutContext(
        MakeCallback(&ConsumerPCC::getdelay, this));
    probeevent = Simulator::ScheduleNow(&ConsumerPCC::PCCProbe, this, 1);
    int i =1;
    for (i=1;i<300;i++)
	Simulator::Schedule(Seconds(i/10.0), &ConsumerPCC::changeprefix, this);
    //Simulator::Schedule(Seconds(2), &ConsumerPCC::changeprefix, this);
    //Simulator::Schedule(Seconds(4), &ConsumerPCC::changeprefix, this);
    //Simulator::Schedule(Seconds(6), &ConsumerPCC::changeprefix, this);
    //Simulator::Schedule(Seconds(8), &ConsumerPCC::changeprefix, this);
    //Simulator::Schedule(Seconds(10), &ConsumerPCC::changeprefix, this);
    //Simulator::Schedule(Seconds(12), &ConsumerPCC::changeprefix, this);
    //Simulator::Schedule(Seconds(14), &ConsumerPCC::changeprefix, this);
    //Simulator::Schedule(Seconds(16), &ConsumerPCC::changeprefix, this);
    //Simulator::Schedule(Seconds(18), &ConsumerPCC::changeprefix, this);
    //Simulator::Schedule(Seconds(20), &ConsumerPCC::changeprefix, this);
    //Simulator::Schedule(Seconds(22), &ConsumerPCC::changeprefix, this);
    //Simulator::Schedule(Seconds(24), &ConsumerPCC::changeprefix, this);
    //Simulator::Schedule(Seconds(26), &ConsumerPCC::changeprefix, this);
    //Simulator::Schedule(Seconds(28), &ConsumerPCC::changeprefix, this);
}

ConsumerPCC::~ConsumerPCC() {}

void ConsumerPCC::changeprefix() {
    NS_LOG_DEBUG(m_interestName);
    if (m_interestName == "/ustc")
        m_interestName = "/ustc1";
    else
        m_interestName = "/ustc";
    NS_LOG_DEBUG(m_interestName);
}

void ConsumerPCC::PCCProbe(uint8_t mode) {
    // NS_LOG_DEBUG("Probe mode: " << (int)mode);
    probestate = mode;
    switch (mode) {
    case 0:
        m_frequency = rate;
        break;
    case 1:
        resettable();
        monitorevent =
            Simulator::Schedule(rtt, &ConsumerPCC::PCCCollect_Calc, this, 0);
        probeevent = Simulator::Schedule(MI, &ConsumerPCC::PCCProbe, this, 2);
        m_frequency = rate * (1 - epsilon);
        break;
    case 2:
        probeevent = Simulator::Schedule(MI, &ConsumerPCC::PCCProbe, this, 0);
        m_frequency = rate * (1 + epsilon);
        break;
    default:
        break;
    }
    NS_LOG_DEBUG("freq: " << m_frequency);
}

void ConsumerPCC::PCCCollect_Calc(uint8_t mode) {
    // NS_LOG_DEBUG("Collect mode: " << (int)mode);
    typedef struct {
        double util;
        double avgrtt;
    } utilret;
    auto calcutil = [this](int i) -> utilret {
        const size_t rttnum = delayMI[i].size();
        assert(rttnum == delayTimestamp[i].size());
        const double lossrate =(double)sqrMI[i].size() / sendnum[i];
        int64_t rttsum = 0;
        double rttgrad = 0;
        double nullval;
        double avgrtt = 0;
        double *x = new double[rttnum];
        double *y = new double[rttnum];
        int index = 0;
     
        while (!delayMI[i].empty()) {
            int64_t rtt = delayMI[i].front().GetNanoSeconds();
            int64_t rttts = delayTimestamp[i].front().GetNanoSeconds();
            x[index] = rttts / 1000000000.0;
            y[index] = rtt / 1000000000.0;
            rttsum += rtt;
            delayMI[i].pop();
            delayTimestamp[i].pop();
            index++;
        }
        gsl_fit_linear(x, 1, y, 1, rttnum, &nullval, &rttgrad, &nullval,
                       &nullval, &nullval, &nullval);
        delete x;
        delete y;
        if (rttgrad < 0.01 || isnan(rttgrad))
            rttgrad = 0;
        NS_LOG_DEBUG("rttsum: " << rttsum << " rttnum: " << rttnum
                                << " rttgrad: " << rttgrad
                                << " lossrate: " << lossrate);

        return {func(rate * (1 + (2 * i - 1) * epsilon), rttgrad, lossrate),
                (double)rttsum / rttnum};
    };

    utilret utils[2];
    double gamma = 0;
    double utildiff = 0;
    double rttdouble = 0;
    double omega = 0;

    switch (mode) {
    case 0:
        monitorevent =
            Simulator::Schedule(MI, &ConsumerPCC::PCCCollect_Calc, this, 1);
        monitorstate = 1;
        break;
    case 1:
        utils[0] = calcutil(0);
        NS_LOG_DEBUG(utils[0].util << ' ' << utils[0].avgrtt);
        monitorevent =
            Simulator::Schedule(MI, &ConsumerPCC::PCCCollect_Calc, this, 2);
        monitorstate = 2;
        break;
    case 2:
        utils[1] = calcutil(1);
        NS_LOG_DEBUG(utils[1].util << ' ' << utils[1].avgrtt);
        utildiff = utils[1].util - utils[0].util;
        rttdouble = utildiff < 0 ? utils[0].avgrtt : utils[1].avgrtt;
        rtt = m_rtt->GetCurrentEstimate();
        NS_LOG_DEBUG("rtt: " << rtt);
        updateMI();
        gamma = m(utildiff) / (2.0 * epsilon * rate);
        omega = (omega0 + omegaboundk * delta);
        if (abs(gamma) > omega * rate) {
            gamma = ((gamma < 0 ? -1 : 1) * omega * rate);
            omegaboundk++;
        } else {
            if (abs(gamma) / rate - omega0 < 0)
                omegaboundk = 0;
            else
                omegaboundk = (int)((abs(gamma) / rate - omega0) / delta) + 1;
        }
        NS_LOG_DEBUG("gamma:" << gamma);
        NS_LOG_DEBUG("omegaboundk:" << omegaboundk);

        rate += gamma;
        NS_LOG_DEBUG("new rate: " << rate);
        if (rate < 100.0)
            rate = 100.0;
        probeevent = Simulator::ScheduleNow(&ConsumerPCC::PCCProbe, this, 1);
        monitorstate = 0;
        break;
    default:
        break;
    }
    m_frequency = rate;
}




double ConsumerPCC::func(double rate, double drtt_dt, double lrate) {
    NS_LOG_DEBUG(rate << " " << drtt_dt << " " << lrate);
    return pow(rate, 0.9) - 900 * rate * drtt_dt - 11.35 * lrate * rate;
}

double ConsumerPCC::m(const double &diff) {
    if (lastdiffpositive ^ (diff > 0)) {
        contk = 1;
        omegaboundk = 0;
    } else {
        contk++;
    }
    lastdiffpositive = (diff > 0);

    NS_LOG_DEBUG("contk:" << contk);

    if (contk <= 3)
        return contk * diff;
    else
        return (2 * contk - 3) * diff;
}

void ConsumerPCC::resettable() {
    for (int i = 0; i < 2; i++) {
        sendnum[i] = 0;
        sqrMI[i].clear();
        while (!delayMI[i].empty()) {
            delayMI[i].pop();
        }
        while (!delayTimestamp[i].empty()) {
            delayTimestamp[i].pop();
        }
    }
}

void ConsumerPCC::updateMI() {
    MI = NanoSeconds(60000000.00);//NanoSeconds(rtt.GetNanoSeconds() * 1.05);
    NS_LOG_DEBUG("MI: " << MI);
}

void ConsumerPCC::getdelay(Ptr<App> app, uint32_t seqno, Time delay,
                           uint32_t retxCount, int32_t hopCount) {
    if (monitorstate <= 0)
        return;
    int index = monitorstate - 1;
    // NS_LOG_DEBUG("try seq " << seqno << " in " << index);
    if (sqrMI[index].find(seqno) != sqrMI[index].end()) {
        sqrMI[index].erase(seqno);
        delayMI[index].emplace(delay);
        delayTimestamp[index].emplace(Simulator::Now());
        // NS_LOG_DEBUG("emplace " << index << " seq " << seqno << " done");
    }
}

void ConsumerPCC::WillSendOutInterest(uint32_t sequenceNumber) {
    // NS_LOG_DEBUG("sendout");
    if (probestate > 0) {
        sqrMI[probestate - 1].emplace(sequenceNumber);
        sendnum[probestate - 1]++;
        // NS_LOG_DEBUG("emplace " << probestate - 1 << " seq " <<
        // sequenceNumber<< " done");
    }
    ScheduleNextPacket();
    ConsumerCbr::WillSendOutInterest(sequenceNumber);
}

void ConsumerPCC::setrate(double ir) {
    rate = ir;
    m_frequency = rate;
}

/*void ConsumerPCC::SetRandomize(const std::string &value)  {
     if (value == "uniform")
     {
           randomSend = CreateObject<UniformRandomVariable>();
           randomSend->SetAttribute("Min", DoubleValue(0.0));
           randomSend->SetAttribute("Max", DoubleValue(2 * 1.0 / frequency));
      }
      else if (value == "exponential")
      {
           randomSend = CreateObject<ExponentialRandomVariable>();
           randomSend->SetAttribute("Mean", DoubleValue(1.0 / frequency));
           randomSend->SetAttribute("Bound", DoubleValue(50.0 / frequency));
      }
      else
           randomSend = 0;
}*/


 void ConsumerPCC::SetRandomPrefix(std::string random_prefix)  {
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



void ConsumerPCC::ScheduleNextPacket()  {
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
	
	ConsumerCbr::ScheduleNextPacket();
}


}  // namespace ndn
}  // namespace ns3
