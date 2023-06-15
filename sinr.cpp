#include "sinr.hpp"

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

NS_LOG_COMPONENT_DEFINE("ndn.sinr");

namespace ns3 {
namespace ndn {

NS_OBJECT_ENSURE_REGISTERED(ConsumerSINR);

TypeId ConsumerSINR::GetTypeId() {
    static TypeId tid =
        TypeId("ns3::ndn::ConsumerSINR")
            .SetGroupName("Ndn")
            .SetParent<ConsumerCbr>()
            .AddConstructor<ConsumerSINR>()
            .AddAttribute("rate", "", DoubleValue(100.0),
                          MakeDoubleAccessor(&ConsumerSINR::setrate),
                          MakeDoubleChecker<double>())
            .AddAttribute("RandomPrefix", "", StringValue(""),
                          MakeStringAccessor(&ConsumerSINR::SetRandomPrefix),
                          MakeStringChecker());

    return tid;
}

ConsumerSINR::ConsumerSINR()
    : ConsumerCbr(),
      rate(0),
      SSflag(1),
      SStime(0),
      fsmax(0.9),
      fsmin(0.5),
      fsold(0.95),
      fsnew(0.95),
      fsstep(0.01),
      Terror(0.2),
      Error_IDG(0),
      IDG_new(0),
      IDG_old(0),
      IIG(0),
      Nidg(20),
      cur_num(20),
      alpha(0.125),
      beta(0.25),
      data_interval(0),
      beforestart(1),
      outtime(0),
      dectime(0),
      outnum(0),
      rttm(0),
      rttnew(0),
      rttold(0),
      varnew(0),
      varold(0),
      RTOnew(0),
      firststartflag(1),
      MI(MilliSeconds(60)){
   
     m_firstInterestDataDelay.ConnectWithoutContext(
        MakeCallback(&ConsumerSINR::getdelay, this));
    // Simulator::Schedule(MI, &ConsumerSINR::ssrate, this );
     
}

ConsumerSINR::~ConsumerSINR() {}

void ConsumerSINR::setrate(double ir) {
    rate = ir;
    if(firststartflag)
    {
        rate = 1000;
        firststartflag=0;
        NS_LOG_DEBUG("first rate: " << rate);
    }

    
    if(rate<200)
        rate=200;
    m_frequency = rate;
    NS_LOG_DEBUG("m_frequency: " << m_frequency);

}

void ConsumerSINR::OnTimeout(uint32_t sequenceNum)
{

    int64_t cur_time = Simulator::Now().GetNanoSeconds();
    outnum+=1;
    NS_LOG_DEBUG("on timeout, sequenceNum: " << sequenceNum);
    if(cur_time - outtime > 60000000)
    {
      
        outtime = Simulator::Now().GetNanoSeconds();
    }
    else
    {
        NS_LOG_DEBUG("on timeout, outnum: " << outnum);
    }
    m_retxSeqs.insert(sequenceNum);
    ScheduleNextPacket();
    
}
void ConsumerSINR::fsupdate()
{

    int64_t cur_time = Simulator::Now().GetNanoSeconds();

    Error_IDG = (IDG_new - IDG_old)/IDG_new;
    NS_LOG_DEBUG("Error_IDG,"<<Error_IDG<<" IDGnew: " << IDG_new << "IDG_old: "<<IDG_old);
    if (Error_IDG <= -Terror && fsnew> fsmin)
    {
        fsnew = fsold - fsstep;
        NS_LOG_DEBUG("rate decrease, fsupdata, fsnew: " << fsnew << "fsold: "<<fsold);
    }
       
    if (outnum>0 && (cur_time - dectime > 60000000))
    {
        fsnew = fsold + 0.5;
        NS_LOG_DEBUG("rate increase, fsupdata, fsnew: " << fsnew << "fsold: "<<fsold);
        dectime = Simulator::Now().GetNanoSeconds();
    }

    

    // if(outnum>1)
    //     fsnew= fsold + fsstep;
    // if(outnum<1)
    //     fsnew= fsold - fsstep;
    
    if (fsnew > 1.5)
        fsnew = 1.5;
    if (fsnew <0.5)
       fsnew = 0.5;   
    
    

}
void ConsumerSINR::ssrate ()
{
     NS_LOG_DEBUG("m_frequency: " << m_frequency);
 	if (SSflag)
	{
        
        NS_LOG_DEBUG("ssrate--enter SSrate: ");
        //rate = m_frequency;
        //m_frequency = rate *  (1+epsilon);
        //rate = m_frequency;
        m_frequency = m_frequency + 10 * SSflag;
        SSflag++;
        NS_LOG_DEBUG("ssrate--m_frequency: " << m_frequency);
            Simulator::Schedule(MI, &ConsumerSINR::ssrate,
                                           this );
	}
}



void ConsumerSINR::resettable() {
    
    data_time.clear();
    cur_num = Nidg;

    fsold = fsnew;
    
    IDG_old = IDG_new;
    data_interval = 0;
    NS_LOG_DEBUG("resettable, fsold " << fsold <<", IDG_old"<<IDG_old<<",datatime.size"<<data_time.size());
}

void ConsumerSINR::getdelay(Ptr<App> app, uint32_t seqno, Time delay,
                                 uint32_t retxCount, int32_t hopCount) {
    
    // NS_LOG_DEBUG("try seq " << seqno << " in " << index);
    rttm = double(delay.GetNanoSeconds());
    rttnew = (1-alpha)*rttold + alpha * rttm;
    varnew = (1-beta)*varold + beta * (std::abs(rttm - rttold));
    RTOnew = rttnew + 4*varnew;
    rttold = rttnew;
    //NS_LOG_DEBUG("RTOnew: "<<RTOnew );
     NS_LOG_DEBUG("m_frequency: " << m_frequency);
       
}

void ConsumerSINR::OnData(shared_ptr<const Data> data) {
    uint64_t congmark = data->getCongestionMark();
    int64_t cur_time = Simulator::Now().GetNanoSeconds();
     ConsumerCbr::OnData(data);
   NS_LOG_DEBUG("Received congestion mark , " << congmark);
    NS_LOG_DEBUG("m_frequency: " << m_frequency);
    if (SSflag)
    {

         NS_LOG_DEBUG("enter slow start: ");
         int64_t cur_time1 = Simulator::Now().GetNanoSeconds();
         if (congmark >2 && (cur_time1 - SStime > 2.5*600000000))
          {
                //NS_LOG_DEBUG("cur_time1: " << cur_time1, << "dectime" << dectime);
                for (int j=0; j< 10; j++)
                     NS_LOG_DEBUG("SS over: ");
               SSflag = 0;
               SStime = Simulator::Now().GetNanoSeconds();
               m_frequency = m_frequency*0.8;
               
          }
         else
         {
            // for (int j=0; j< 10; j++)
              //   NS_LOG_DEBUG("enter SS: ");
             
	     //
             m_frequency=m_frequency+5;
             
             NS_LOG_DEBUG("mark: " << congmark);
             //NS_LOG_DEBUG("SSflag: "<<SSflag);
         }
    }
   
    if(!SSflag)
   {

        if(cur_num >=0)
        {
            int64_t cur_time2 = Simulator::Now().GetNanoSeconds();
            data_time.push_back(cur_time2);
            cur_num--;
        }
        else
        {
            for(int i =0; i<Nidg; i++)
            {
                int64_t a = data_time[i+1]-data_time[i];
                NS_LOG_DEBUG("data_interval: " << a <<"index i"<<i);
                data_interval += a;
            }
            if(beforestart >0 )
            {
                IDG_new = double(data_interval)/Nidg;
                NS_LOG_DEBUG("IDG_new: " << IDG_new);
                beforestart--;
            }
            else
            {

                double b =double(m_rtt->GetCurrentEstimate().GetNanoSeconds())*outnum;
                NS_LOG_DEBUG("timeout delay: " << b);
                IDG_new = double(data_interval)/Nidg;
                
                NS_LOG_DEBUG("IDG_old: " << IDG_old<<"IDG_new: " << IDG_new);
                // if(outnum==0)
                // {
                //     if(IDG_new<=IDG_old )
                //     {
                //         if(fsnew>1.2)
                //         fsnew=1.2;
                //         fsnew-=0.01;
                //         NS_LOG_DEBUG("fs decrease: ");
                //     }
                  
                //     if(fsnew<0.5)
                //         fsnew=0.5;
                //     NS_LOG_DEBUG("fsold: " << fsold<<"fsnew: " << fsnew);
                // }
                // else
                // {
                //     if(cur_time-dectime  > 6000000)
                //     {
                //         if(fsnew<1)
                //             fsnew=1;
                //          fsnew+=0.2;
                //          NS_LOG_DEBUG("fs increase: ");
                //          if(fsnew>1.5)
                //              fsnew=1.5;
                //          NS_LOG_DEBUG("fsold: " << fsold<<"fsnew: " << fsnew);
                //     }
                //     dectime = Simulator::Now().GetNanoSeconds();

                 
                // } 
                if(IDG_new >IIG)
                {
                    IIG=1.03*IDG_new;
                    NS_LOG_DEBUG("IDG_new > IIG, congestion occurs, IDG_new: "<< IDG_new<<"IIG update, IIG: " << IIG);   
                }
                else if(IDG_new<=IIG)
                {
                    IIG = 0.95*IDG_new;
                     NS_LOG_DEBUG("IDG_new ==IIG, path less utilized, increase rate, IDG_new: "<< IDG_new<<"IIG update, IIG: " << IIG);  
                }
                //NS_LOG_DEBUG("fsnew: " << fsnew<<" IDG_new: " << IDG_new<<" IIG:"<<fsnew*IDG_new);               
                //IIG =fsnew*IDG_new;    
                 NS_LOG_DEBUG("IDG: "<< IDG_new<<"IIG: " << IIG);            
                NS_LOG_DEBUG("IIG: " << IIG <<" RTOnew: "<<RTOnew<<" RTOsystem:"<<m_rtt->RetransmitTimeout().GetNanoSeconds());
                setrate((1.0 * 1000000000)/IIG);
                outnum = 0;
            
            }
            resettable();
            
        }
        
        
    }
    
     
	    
}
    
    //auto da = datainfo[index].find(seqno);

void ConsumerSINR::SetRandomPrefix(std::string random_prefix) {
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

void ConsumerSINR::ScheduleNextPacket() {
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
/*
void ConsumerSINR::changeprefix() {
    NS_LOG_DEBUG(m_interestName);
    if (m_interestName == "/ustc")
        m_interestName = "/ustc1";
    else
        m_interestName = "/ustc";
    NS_LOG_DEBUG(m_interestName);
}*/

}  // namespace ndn
}  // namespace ns3
