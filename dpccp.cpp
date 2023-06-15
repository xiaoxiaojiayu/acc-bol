#include "dpccp.hpp"

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

namespace ns3 {
namespace ndn {

NS_LOG_COMPONENT_DEFINE("ndn.DPCCP");

enum class ConsumerDPCCP::States { READY, RTTPROBE, BF, FC, QD };
enum class ConsumerDPCCP::BFStates { PLUS, NOCHANGE, MINUS, DECISION };
enum class ConsumerDPCCP::FCStates { S0, S1, S2, DECISION };

NS_OBJECT_ENSURE_REGISTERED(ConsumerDPCCP);

std::ostream& operator<<(std::ostream& os, ConsumerDPCCP::States s) {
    static const std::string str[5] = {"READY", "RTTPROBE", "BF", "FC", "QD"};
    os << str[(int)s];
    return os;
}

std::ostream& operator<<(std::ostream& os, ConsumerDPCCP::BFStates s) {
    static const std::string str[4] = {"PLUS", "NOCHANGE", "MINUS", "DECISION"};
    os << str[(int)s];
    return os;
}

std::ostream& operator<<(std::ostream& os, ConsumerDPCCP::FCStates s) {
    static const std::string str[4] = {"S0", "S1", "S2", "DECISION"};
    os << str[(int)s];
    return os;
}

ConsumerDPCCP::ConsumerDPCCP()
    : curstate(States::READY),
      beta(0.05),
      zeta(0.8),
      psi(3),
      alpha_r(20),
      rho(1),
      eta(30),
      settime0(0),
      settime1(0),
      srcnum(3),
      SSflag(0),
      rttupdateflag(0) {
    m_lastRetransmittedInterestDataDelay.ConnectWithoutContext(
        MakeCallback(&ConsumerDPCCP::getdelay, this));
    NS_LOG_DEBUG("Enter state: " << curstate);
    baseRtt[0] = MilliSeconds(60);
    baseRtt[1] = MilliSeconds(60);
    baseRtt[2] = MilliSeconds(60);
    baseRtt0_step=NanoSeconds(100000000);
    baseRtt1_step=NanoSeconds(100000000);
    maxrtt = MilliSeconds(60);
    
};

ConsumerDPCCP::~ConsumerDPCCP() {}

void ConsumerDPCCP::StartApplication(){
    Consumer::StartApplication();
    nextevent = ns3::Simulator::ScheduleNow(&ConsumerDPCCP::ProbeRtt, this);
}

TypeId ConsumerDPCCP::GetTypeId() {
    static TypeId tid = TypeId("ns3::ndn::ConsumerDPCCP")
                            .SetGroupName("Ndn")
                            .SetParent<ConsumerCbr>()
                            .AddConstructor<ConsumerDPCCP>()
    			    .AddAttribute(
                "RandomPrefix", "", StringValue(""),
                MakeStringAccessor(&ConsumerDPCCP::SetRandomPrefix),
                MakeStringChecker());	
    return tid;
}

void ConsumerDPCCP::ProbeRtt() {
    curstate = States::RTTPROBE;
    NS_LOG_DEBUG("Enter ProneRTT: ");
    int64_t cur_time1 = Simulator::Now().GetNanoSeconds();
    if(cur_time1 > 1000000000)
        rate = m_frequency = 100;
    else
        rate = m_frequency;
    //baseRtt0 = MilliSeconds(20);
    //baseRtt1 = MilliSeconds(60);
     NS_LOG_DEBUG("baseRtt[0] :" << baseRtt[0]);
    
    MI = NanoSeconds(maxrtt.GetNanoSeconds());
    //nextevent.Cancel();
    NS_LOG_DEBUG("MI :" << MI);
    nextevent = ns3::Simulator::ScheduleNow(&ConsumerDPCCP::BFEntry, this);//ns3::MilliSeconds(10),
}

void ConsumerDPCCP::BFEntry() {
    curstate = States::BF;
    rate = m_frequency;
    NS_LOG_DEBUG("Enter state: " << curstate);
    nextevent = ns3::Simulator::ScheduleNow(&ConsumerDPCCP::BFDoing, this,
                                            BFStates::PLUS);
}

void ConsumerDPCCP::BFDoing(BFStates bfs) {
    NS_LOG_DEBUG("curstate: " << curstate << "::" << bfs);
    switch (bfs) {
    case BFStates::PLUS:
        m_frequency = (1 + beta) * rate;
        nextevent = ns3::Simulator::Schedule(MI, &ConsumerDPCCP::BFDoing, this,
                                             BFStates::NOCHANGE);
        break;
    case BFStates::NOCHANGE:
        m_frequency = rate;
        nextevent = ns3::Simulator::Schedule(MI, &ConsumerDPCCP::BFDoing, this,
                                             BFStates::MINUS);
        break;
    case BFStates::MINUS:
        m_frequency = (1 - beta) * rate;
        nextevent = ns3::Simulator::Schedule(MI, &ConsumerDPCCP::BFDoing, this,
                                             BFStates::DECISION);
        break;
    case BFStates::DECISION:
        auto ret = calc_backlogged_by_datainfomap();
        double backlogged = ret.backlogged;

        if (backlogged <= psi) {
            rate = (1 + beta) * rate;
            NS_LOG_DEBUG("still in BF, rate=" << rate);
            nextevent = ns3::Simulator::ScheduleNow(&ConsumerDPCCP::BFDoing,
                                                    this, BFStates::PLUS);
        } else {
            rate = ret.receive_rate;
            fc_delta_r[2] = backlogged;
            NS_LOG_DEBUG("switch to FC, rate=" << rate);
            nextevent =
                ns3::Simulator::ScheduleNow(&ConsumerDPCCP::FCEntry, this);
        }
        break;
    }
}

   
void ConsumerDPCCP::FCEntry() {
    curstate = States::FC;
    NS_LOG_DEBUG("Enter state: " << curstate);
    rho = 1;
    fc_delta_r[0] = fc_delta_r[1] = fc_delta_r[2];
    nextevent = ns3::Simulator::ScheduleNow(&ConsumerDPCCP::FCDoing, this,
                                            FCStates::S0);
}

void ConsumerDPCCP::FCDoing(FCStates fcs) {
    NS_LOG_DEBUG("curstate: " << curstate << "::" << fcs);
    switch (fcs) {
    case FCStates::S0:
        if (fc_delta_r[2] < alpha_r)
            rate = rate + rho;
        if (fc_delta_r[2] > alpha_r)
            rate = rate - rho;
        m_frequency = rate;
        nextevent = ns3::Simulator::Schedule(MI, &ConsumerDPCCP::FCDoing, this,
                                             FCStates::S1);
        break;
    case FCStates::S1:
        fc_delta_r[0] = calc_backlogged_by_datainfomap().backlogged;
        nextevent = ns3::Simulator::Schedule(MI, &ConsumerDPCCP::FCDoing, this,
                                             FCStates::S2);
        break;
    case FCStates::S2:
        fc_delta_r[1] = calc_backlogged_by_datainfomap().backlogged;
        nextevent = ns3::Simulator::Schedule(MI, &ConsumerDPCCP::FCDoing, this,
                                             FCStates::DECISION);
        break;
    case FCStates::DECISION:
        fc_delta_r[2] = calc_backlogged_by_datainfomap().backlogged;
        double backlogged = fc_delta_r[2];
        if (backlogged <= psi) {
            nextevent =
                ns3::Simulator::ScheduleNow(&ConsumerDPCCP::BFEntry, this);
        } else if (backlogged > eta) {
            nextevent =
                ns3::Simulator::ScheduleNow(&ConsumerDPCCP::QDEntry, this);
        } else {
            double abs_fc_delta_r_rate =
                std::abs((fc_delta_r[2] - fc_delta_r[0]) / 3);
            if (abs_fc_delta_r_rate < 1)
                rho++;
            if (abs_fc_delta_r_rate > 1)
                rho--;
            if (fc_delta_r[2] >= alpha_r && alpha_r >= fc_delta_r[1])
                rho = 1;
            if (fc_delta_r[2] <= alpha_r && alpha_r <= fc_delta_r[1])
                rho = 1;
            nextevent = ns3::Simulator::ScheduleNow(&ConsumerDPCCP::FCDoing,
                                                    this, FCStates::S0);
        }
        break;
    }
}

void ConsumerDPCCP::QDEntry() {
    curstate = States::QD;
    NS_LOG_DEBUG("Enter state: " << curstate);
    m_frequency = MIN_RATE;
    nextevent = ns3::Simulator::Schedule(MI / 2, &ConsumerDPCCP::QDDoing, this);
}

void ConsumerDPCCP::QDDoing() {
    double backlogged = calc_backlogged_by_datainfomap().backlogged;
    if (backlogged != 0.0)
        nextevent =
            ns3::Simulator::Schedule(MI / 2, &ConsumerDPCCP::QDDoing, this);
    else {
        m_frequency = zeta * rate;
        rate = m_frequency;
        nextevent = ns3::Simulator::ScheduleNow(&ConsumerDPCCP::BFEntry, this);
    }
}

int64_t ConsumerDPCCP::getSrcofData(shared_ptr<const Data> data) {
    auto blk =
        data->getMetaInfo().findAppMetaInfo(::ndn::tlv::AppPrivateBlock1);
    if (blk == nullptr)
        return -1;
    else
        return ::ndn::readNonNegativeInteger(*blk);
}

void ConsumerDPCCP::OnData(shared_ptr<const Data> data) {
    uint32_t seqNum = data->getName().get(-1).toSequenceNumber();
    int64_t respondid = getSrcofData(data);
    //NS_LOG_DEBUG("respondid " << respondid );

    datainfomap.emplace(seqNum, datainfo{seqNum, respondid, 0, 0});

    ConsumerCbr::OnData(data);
}

void ConsumerDPCCP::getdelay(Ptr<App> app, uint32_t seqno, Time delay,
                             int32_t hopCount) {
    // NS_LOG_DEBUG("try seq " << seqno << " in " << index);
    auto it = datainfomap.find(seqno);
    if (it != datainfomap.end()) {
        it->second.rtt = delay.GetNanoSeconds();
        it->second.timestamp = Simulator::Now().GetNanoSeconds();
    }
}

void ConsumerDPCCP::ScheduleNextPacket() {
    if (randprefix) {
        static std::random_device rdev;
        static std::mt19937 reng(rdev());
        static std::uniform_real_distribution<> u(0, intervalp.back());
        double rnum = u(reng);
        int i;
        for (i = 0; i < intervalp.size(); i++)
            if (rnum <= intervalp[i])
                break;
        this->SetAttribute("Prefix",
                           StringValue(intervalmap.find(intervalp[i])->second));
    }
    ConsumerCbr::ScheduleNextPacket();
}

void ConsumerDPCCP::SetRandomPrefix(std::string random_prefix) {
    if (!random_prefix.empty()) {
        double acml = 0;
        std::size_t bpos = 0, epos = 0;
        do {
            epos = random_prefix.find(',', bpos);
            std::string kvpair = random_prefix.substr(bpos, epos - bpos);
            bpos = epos + 1;
            acml += std::stod(kvpair.substr(kvpair.find('=') + 1));
            intervalp.push_back(acml);
            intervalmap.emplace(
                std::make_pair(acml, kvpair.substr(0, kvpair.find('='))));
        } while (bpos != std::string::npos + 1);
        randprefix = true;
    } else
        randprefix = false;
}

ConsumerDPCCP::bcrettype ConsumerDPCCP::calc_backlogged_by_datainfomap() {
    int64_t rttsum[srcnum] = {0, 0,0}, rttnum[srcnum] = {0, 0,0};
    
    maxrtt = NanoSeconds(100000);
    for (auto& i : datainfomap) {
        int index = 0;
        NS_LOG_DEBUG("i.second.rtt: " << i.second.rtt);
        NS_LOG_DEBUG("maxrtt: " << maxrtt);

         if( NanoSeconds(i.second.rtt) > maxrtt && i.second.rtt < 200000000)
        { 
            maxrtt = NanoSeconds(i.second.rtt);
        }
       
        if (i.second.respondid==0)
        {
            index= 0;
            //NS_LOG_DEBUG("i.second.rtt: " << i.second.rtt);
            if(NanoSeconds(i.second.rtt) < baseRtt[index])
            {
                
                
                baseRtt[index] = NanoSeconds(i.second.rtt);                 

                 NS_LOG_DEBUG("find new baseRtt1: " << baseRtt[index]);
            
             }
        }
            
        else if(i.second.respondid==1)
        {
            index =1;
            if(NanoSeconds(i.second.rtt) < baseRtt[index])
            {
                 baseRtt[index] = NanoSeconds(i.second.rtt);
                 
                 NS_LOG_DEBUG("find new baseRtt0: " << baseRtt[index]);
            }

        }
        else if(i.second.respondid==2 )
        {
            index =2;
            if(NanoSeconds(i.second.rtt) < baseRtt[index])
            {
                 baseRtt[index] = NanoSeconds(i.second.rtt);
                 
                 NS_LOG_DEBUG("find new baseRtt0: " << baseRtt[index]);
            }

        }
        rttupdateflag =1;
        // int64_t cur_time = Simulator::Now().GetNanoSeconds();
        // if((cur_time - settime0 > 5*60000000))
        // {
        //         baseRtt0 = baseRtt0_step;
        //         baseRtt0_step=NanoSeconds(100000000);
        //         settime0 = Simulator::Now().GetNanoSeconds();
        //         NS_LOG_DEBUG("set basertt0: " << baseRtt1);
        // }
            
        // else
        //     baseRtt0 = std::min(baseRtt0_step,baseRtt0);
        // if((cur_time - settime1 > 5*60000000))
        // {
        //     baseRtt1 = baseRtt1_step;
        //     baseRtt1_step=NanoSeconds(200000000);
        //     settime1 = Simulator::Now().GetNanoSeconds();
        //     NS_LOG_DEBUG("set basertt1: " << baseRtt1);
        // }
            
        // else
        // {
        //     baseRtt1 = std::min(baseRtt1_step, baseRtt1);
        //     //baseRtt1_step=NanoSeconds(100000000);
        // }
            
        
            
        rttsum[index] += i.second.rtt;
        rttnum[index]++;

         
        
    }
    NS_LOG_DEBUG("baseRtt0,"<<baseRtt[0]);
    NS_LOG_DEBUG("baseRtt1,"<<baseRtt[1]);
    NS_LOG_DEBUG("baseRtt2,"<<baseRtt[2]);
    NS_LOG_DEBUG("rttsum{src0,src1,src2} = {" << rttsum[0] << "," << rttsum[1] << "," << rttsum[2]
                                         << "}");
    NS_LOG_DEBUG("rttnum{src0,src1, src2} = {" << rttnum[0] << "," << rttnum[1] << "," << rttnum[2]
                                         << "}");

    int64_t div[srcnum] = {rttnum[0], rttnum[1], rttnum[2]};
    int64_t number_of_MI = 3;

    for (int i = 0; i < srcnum; i++)
        if (div[i] == 0)
            div[i]++;
    if (curstate == States::FC) {
        number_of_MI = 1;
    }

    double backlogged = 0;
    for(int i=0; i<srcnum; i++)
    {
        backlogged+=(double)rttnum[i] / number_of_MI / MI.GetNanoSeconds() *
            ((double)(rttsum[i] - baseRtt[i].GetNanoSeconds() * rttnum[i]) /
             div[i]);

    NS_LOG_DEBUG("rttnum[i],"<<rttnum[i]<<"number_of_MI:"<<number_of_MI<<"{MI.GetNanoSeconds()}:"<<MI.GetNanoSeconds()
                <<"{rttsum[i] - baseRtt1.GetNanoSeconds() * rttnum[i]}:"<<(double)(rttsum[i] - baseRtt1.GetNanoSeconds() * rttnum[i])<<"div[i]:"<<div[i]);
    }
    double receive_rate =0;
    for(int i=0; i<srcnum; i++)
    {
      receive_rate+= rttnum[i]*
                          (1000000000.0 / number_of_MI) / MI.GetNanoSeconds();
    }
    datainfomap.clear();
    NS_LOG_DEBUG("{backlogged,receive_rate} = {" << backlogged << ","
                                                 << receive_rate << "}");
    return {backlogged, receive_rate};
}

}  // namespace ndn
}  // namespace ns3
