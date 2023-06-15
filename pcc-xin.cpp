#include "pcc-xin.hpp"

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

NS_LOG_COMPONENT_DEFINE("ndn.PCCmy");

namespace ns3 {
namespace ndn {

NS_OBJECT_ENSURE_REGISTERED(ConsumerPCCMY);

TypeId ConsumerPCCMY::GetTypeId() {
    static TypeId tid =
        TypeId("ns3::ndn::ConsumerPCCMY")
            .SetGroupName("Ndn")
            .SetParent<ConsumerCbr>()
            .AddConstructor<ConsumerPCCMY>()
            .AddAttribute("rate", "", DoubleValue(3000.0),
                          MakeDoubleAccessor(&ConsumerPCCMY::setrate),
                          MakeDoubleChecker<double>())
            .AddAttribute("epsilon", "", DoubleValue(0.05),
                          MakeDoubleAccessor(&ConsumerPCCMY::epsilon),
                          MakeDoubleChecker<double>())
            .AddAttribute("RandomPrefix", "", StringValue(""),
                          MakeStringAccessor(&ConsumerPCCMY::SetRandomPrefix),
                          MakeStringChecker());

    return tid;
}

ConsumerPCCMY::ConsumerPCCMY()
    : ConsumerCbr(),
      rate(0),
      rtt(MilliSeconds(70)),
      MI(MilliSeconds(70)),
      minrtt(MilliSeconds(20)),
      maxrtt(MilliSeconds(70)),
      rttdiff(MilliSeconds(0)),
      epsilon(0.05),
      theta0(1),
      omega0(0.05),
      delta(0.1),
      omegaboundk(0),
      probestate(0),
      monitorstate(0),
      lastdiffpositive(false),
      contk(1),
      SSflag(1),
      SStime(0),
      dectime(0),
      fincreaseflag(0),
      insmarkrate(0),
      heavycong(0),
      decflag(1),
      receivenum(0),
      mrate0(0),
      srcnum(3),
      timeoutnum(0) {
    m_firstInterestDataDelay.ConnectWithoutContext(
        MakeCallback(&ConsumerPCCMY::getdelay, this));
    
    maxrtt = NanoSeconds(70000000);
    minrtt = NanoSeconds(20000000);

    //   Simulator::Schedule(MI, &ConsumerPCCMY::ssrate, this );
     //for (int i = 1; i < 60; i++)
    // Simulator::Schedule(Seconds(i / 2.0), &ConsumerPCCMY::changeprefix, this);
}

ConsumerPCCMY::~ConsumerPCCMY() {}

void ConsumerPCCMY::StartApplication(){
    Consumer::StartApplication();
    Simulator::Schedule(MI, &ConsumerPCCMY::ssrate, this );
}

void ConsumerPCCMY::changeprefix() {
    NS_LOG_DEBUG(m_interestName);
    if (m_interestName == "/ustc")
        m_interestName = "/ustc1";
    else
        m_interestName = "/ustc";
    NS_LOG_DEBUG(m_interestName);
}

void ConsumerPCCMY::resettable() {
    for (int i = 0; i < 3; i++) {
        sendnum[i] = 0;
        sqrMI[i].clear();
        datainfo[i].clear();
    }
    lossnum = 0;
    marknum = 0;
    marksum = 0;
    insmarkrate = 0;
    receivenum = 0;
    mrate0 = 0;
    timeoutnum = 0;
    //lossrate = 0;
}

void ConsumerPCCMY::updateMI() {
    //MI = NanoSeconds(rtt.GetNanoSeconds() * 1.05);
    MI=NanoSeconds(maxrtt.GetNanoSeconds());
    NS_LOG_DEBUG("MI: " << MI);
}

void ConsumerPCCMY::setrate(double ir) {
    rate = ir;
    m_frequency = rate;

}

void ConsumerPCCMY::PCCProbe(uint8_t mode) {
    // NS_LOG_DEBUG("Probe mode: " << (int)mode);
    probestate = mode;
    switch (mode) {
    case ENTER_PROBE_RATE:
        m_frequency = rate;
        break;
    case ENTER_PROBE_RATE_PLUS:
        resettable();
        monitorevent = Simulator::Schedule(rtt, &ConsumerPCCMY::PCCCollect_Calc,
                                           this, ENTER_MONITOR_RATE_PLUS);
        probeevent = Simulator::Schedule(MI, &ConsumerPCCMY::PCCProbe, this,
                                         ENTER_PROBE_RATE_MINUS);
        m_frequency = rate * (1 - epsilon);
        break;
    case ENTER_PROBE_RATE_MINUS:
        probeevent = Simulator::Schedule(MI, &ConsumerPCCMY::PCCProbe, this,
                                         ENTER_PROBE_RATE);
        m_frequency = rate * (1 + epsilon);
        break;
    default:
        break;
    }
    NS_LOG_DEBUG("freq: " << m_frequency);
}

double ConsumerPCCMY::funcA(double rate, double drtt_dt1, double weight1,
                           double drtt_dt2, double weight2, double markrate) {
    NS_LOG_DEBUG(rate << " " << drtt_dt1 << " " << drtt_dt2 << " " << markrate);
    double w =0;
    double maxcong = 0;
    if (drtt_dt1 > drtt_dt2)
    {
	maxcong = drtt_dt1;
        w = weight1;
        return pow(rate, 0.9) - 900 * rate * weight1 * drtt_dt1 -11.35 * markrate * rate * weight1;
    }
     else
    {
	maxcong = drtt_dt2;
        w = weight2;
        return pow(rate, 0.9) - 900 * rate * weight2 * drtt_dt2 -11.35 * markrate * rate * weight2;
    }

    //return pow(rate, 0.9) - 900 * rate * weight2 * drtt_dt2 -
    //       11.35 * markrate * rate * weight2;
}

double ConsumerPCCMY::func1(double rate, double markrate, double lossrate) {
    // NS_LOG_DEBUG(rate << " " << drtt_dt1 << " " << drtt_dt2 << " " <<
    // markrate);
    return 1.1*pow(rate, 0.9) - 900 * markrate * rate - 11.35 * lossrate * rate;
}

double ConsumerPCCMY::funcB(double rate, double drtt_dt1,   double weight1, double drtt_dt2,  double weight2, double markrate) {
    NS_LOG_DEBUG(rate << " " << drtt_dt1 << " " << drtt_dt2 << " " << markrate);
    //return pow(rate, 0.9) - 100 * rate * (weight1 * drtt_dt1 + weight2 * drtt_dt2) - 5.35 * markrate * rate;// ;//
    //return pow(rate, 0.9) - 900 * rate * (weight1 * drtt_dt1 + weight2 * drtt_dt2) - 11.35 * markrate * rate;
    return pow(rate, 0.9) -900* rate*drtt_dt1 -11.35 * rate*markrate; 
           // 10 * rate * (weight1 * drtt_dt1 + weight2 * drtt_dt2) -
}

double ConsumerPCCMY::funcC(double rate, double drtt_dt1,   double weight1, double drtt_dt2,  double weight2, double lrate) {
    NS_LOG_DEBUG(rate << " " << drtt_dt1 << " " << drtt_dt2 << " " << lrate);
    return pow(rate, 0.9) - 900 * rate * (weight1 * drtt_dt1 + weight2 * drtt_dt2) - 11.35 * lrate * rate;
}

void ConsumerPCCMY::ssrate ()
{
 	if (SSflag)
	{
        
        NS_LOG_DEBUG("ssrate--enter SSrate: ");
        //rate = m_frequency;
        //m_frequency = rate *  (1+epsilon);
        //rate = m_frequency;
        m_frequency = m_frequency + 20 * SSflag;
        SSflag++;
        NS_LOG_DEBUG("ssrate--m_frequency: " << m_frequency);
            Simulator::Schedule(MI, &ConsumerPCCMY::ssrate,
                                           this );
	}
}

double ConsumerPCCMY::m(const double &diff) {
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
        return (2.5 * contk - 3) * diff;
}

utilret ConsumerPCCMY::calcutil(int i) {
    NS_LOG_DEBUG("entercalcutil: ");
    i++;
    const size_t rttnum = datainfo[i].size();
    // double lossrate=0;
    // lossrate = (double)sqrMI[i].size() / (sendnum[i] + 1);
    double lossrate=double(lossnum) / rttnum;  //
    // if (lossrate < 0.009)
    //     lossrate = 0;
    const double markrate = double(marknum) / receivenum;  //(sendnum[i]+1);
    

    NS_LOG_DEBUG("rttnum: " << rttnum);
    NS_LOG_DEBUG("sendnum[i]: " << sendnum[i]);
    NS_LOG_DEBUG("(double)sqrMI[i].size(): " << (double)sqrMI[i].size());
    //NS_LOG_DEBUG("lossrate: " << lossrate);
    NS_LOG_DEBUG("marknum: " << marknum << "receivenum: " << receivenum);
    NS_LOG_DEBUG("markrate: " << markrate);

    double nullval;
    const size_t srcnum = 3;
    double *avgrtt= new double[srcnum];
    double *weight= new double[srcnum];
    double *mrate= new double[srcnum];
    int *recenum= new int[srcnum];
    int *mksum= new int[srcnum];
    int *mknum= new int[srcnum];
    double *smrtt= new double[srcnum];
    double *rttgrade = new double[srcnum];
    for (int src = 0; src < srcnum; src++)
    {
        avgrtt[src]=0;
        weight[src]=0;
        mrate[src]=0;
        recenum[src]=0;
        mksum[src]=0;
        mknum[src]=0;
        smrtt[src]=0;
        rttgrade[src]=0;
    }
    int64_t rttsum = 0;

    int64_t rttsum1 = 0;
    double rttgrad1 = 0;
    double avgrtt1 = 0;
    int rttnum1 = 0;
    double weight1 = 0;
    double markrate1 = 0;
    int64_t rttsum2 = 0;
    double rttgrad2 = 0;
    double avgrtt2 = 0;
    int rttnum2 = 0;
    double weight2 = 0;
    double avgmark1 =0;
    double avgmark2 =0;
    int rttnum3 = 0;
    double rttgrad3 = 0;

    double *x1 = new double[rttnum];
    double *y1 = new double[rttnum];
    double *x2 = new double[rttnum];
    double *y2 = new double[rttnum];
    double *x3 = new double[rttnum];
    double *y3 = new double[rttnum];

    auto it = datainfo[i].cbegin();

    while (it != datainfo[i].cend()) {

        if(NanoSeconds(it->second.rtt) > maxrtt && it->second.rtt < 200000000)
        {
            
            maxrtt = NanoSeconds(it->second.rtt);
            
        }   
        if(NanoSeconds(it->second.rtt) < minrtt)
        {
            minrtt = NanoSeconds(it->second.rtt);
        }

        int cursrc = it->second.src;
        smrtt[cursrc] += it->second.rtt / 1000000.0;
        //NS_LOG_DEBUG(" it->second.rtt" << it->second.rtt / 1000000.0);
        recenum[cursrc]+=1;
       

        if(it->second.congmark >30)
        {
            mknum[cursrc]+=1;

            mksum[cursrc]+=it->second.congmark ;
           // mrate0=1;
        }

        // NS_LOG_DEBUG("enter while circle " );
        int64_t rtt = it->second.rtt;
        int64_t ts = it->second.timestamp;
        if (cursrc==0) {
            // NS_LOG_DEBUG("The source is r0: " );
            x1[rttnum1] = ts / 1000000000.0;
            y1[rttnum1] = rtt / 1000000000.0;
            //y1[rttnum1] = it->second.congmark;
            rttsum1 += rtt;
            rttnum1++;
            rttsum += rtt;
        } else if (cursrc==1) {
            // NS_LOG_DEBUG("The source is p0: " );
            x2[rttnum2] = ts / 1000000000.0;
            y2[rttnum2] = rtt / 1000000000.0;
            //y2[rttnum2] = it->second.congmark;
            rttsum2 += rtt;
            rttnum2++;
            rttsum += rtt;
        }else if(cursrc==2)
        {
            x3[rttnum3] = ts / 1000000000.0;
            y3[rttnum3] = rtt / 1000000000.0;
            rttnum3++;
        }

        it++;
    }
    for (int src = 0; src < srcnum; src++)
    {
        avgrtt[src]=smrtt[src]/(1+recenum[src]);
        weight[src]=double(recenum[src])/(1+rttnum);
        mrate[src]=double(mknum[src])/(1+recenum[src]);
    }

    
    

    double temp=0;
    for(int j=0; j<rttnum1; j++)
    {
        temp+=y1[j];
        //NS_LOG_DEBUG("y1[j]: " << y1[j]);
    }
    avgmark1 = double(temp)/(32.0*(rttnum1+1));
    NS_LOG_DEBUG("avgmark1: " << avgmark1);
    temp=0;
     for(int j=0; j<rttnum2; j++)
    {
        temp+=y2[j];
        //NS_LOG_DEBUG("y1[j]: " << y2[j]);
    }

    avgmark2 = double(temp)/(32.0* (rttnum2+1));
    NS_LOG_DEBUG("avgmark2: " << avgmark2);
    

    gsl_fit_linear(x1, 1, y1, 1, rttnum1, &nullval, &rttgrade[0], &nullval,
                   &nullval, &nullval, &nullval);
    gsl_fit_linear(x2, 1, y2, 1, rttnum2, &nullval, &rttgrade[1], &nullval,
                   &nullval, &nullval, &nullval);

    gsl_fit_linear(x3, 1, y3, 1, rttnum3, &nullval, &rttgrade[2], &nullval,
                   &nullval, &nullval, &nullval);


    double avgcongestion = 0;
    for (int src = 0; src < srcnum; src++)
    {
        avgcongestion+=rttgrade[src]*weight[src];
    } 
   mrate0 = double(mknum[0]+mknum[1]+mknum[2])/(1+recenum[0]+recenum[1]+recenum[2]);
     
    delete[] x1;
    delete[] y1;
    delete[] x2;
    delete[] y2;

    if (abs(rttgrad1) < 0.01 || isnan(rttgrad1))
        rttgrad1 = 0;
    if (abs(rttgrad2) < 0.01 || isnan(rttgrad2))
        rttgrad2 = 0;
    weight1 = double(rttnum1) / (rttnum+1);
    weight2 = double(rttnum2) / (rttnum+1);
    // NS_LOG_DEBUG("rttsum1: " << rttsum1 << " rttnum1: " << rttnum1
    //                          << " weight1: " << weight1 << " markgrad1: "
    //                          << rttgrad1 << " markrate: " << markrate);
    // NS_LOG_DEBUG("rttsum2: " << rttsum2 << " rttnum2: " << rttnum2
    //                          << " weight2: " << weight2 << " markgrad2: "
    //                          << rttgrad2 << " markrate: " << markrate);
    // NS_LOG_DEBUG("avgmark1: " << avgmark1 << " recenum[0]: " << recenum[0]
    //                          << " weight1: " << weight[0] << " mknum[0]: "
    //                          << mknum[0] <<" mkrate[0]: " << mrate[0]);
    // NS_LOG_DEBUG("avgmark2: " << avgmark2 << " recenum[1]: " << recenum[1]
    //                          << " weight2: " << weight[1]<< " mknum: "
    //                          << mknum[1]<<" mkrate[1]: " << mrate[1]  );

    for(int i=0;i<srcnum ; i++)
    {
    NS_LOG_DEBUG("rttgrade: " << rttgrade[i] 
                             << " weight: " << weight[i] 
                             <<"recenum: "<<recenum[i]);
    }

    i--;
    //mrate0 = std::max(mrate[0],mrate[1]);
    NS_LOG_DEBUG("mrate0: " << mrate0);
    /*if(mrate0 < 0.001)
    {
        if(randomloss!=0 && randomloss>lossrate)
            randomloss=lossrate;
        if(randomloss==0)
            randomloss=lossrate;
        
    }*/
    //lossrate=0;
    //NS_LOG_DEBUG("randomloss: " << randomloss);
    // return {funcA(rate * (1 + (2 * i - 1) * epsilon), rttgrad1, weight1,
    //              rttgrad2, weight2, markrate),
    //         (double)rttsum / rttnum};
    return {funcB(rate * (1 + (1 - 2 * i) * epsilon), avgcongestion, weight1,
                  avgmark2, weight2, mrate0)};//mrate[0]
    //return {funcC(rate * (1 + (2 * i - 1) * epsilon), rttgrad1, weight1,
    //             rttgrad2, weight2, lossrate),
    //        (double)rttsum / rttnum};

    // return {func1(rate * (1 + (2 * i - 1) * epsilon), markrate, lossrate),
    //         (double)rttsum / rttnum};
   
}


void ConsumerPCCMY::PCCCollect_Calc(uint8_t mode) {
    NS_LOG_DEBUG("Collect_calculate mode: " << (int)mode);

    utilret utils[2];
    double gamma = 0;
    double utildiff = 0;
    double rttdouble = 0;
    double omega = 0;

    switch (mode) {
    case ENTER_MONITOR_RATE_PLUS:
        monitorevent = Simulator::Schedule(MI, &ConsumerPCCMY::PCCCollect_Calc,
                                           this, ENTER_MONITOR_RATE_MINUS);
        monitorstate = 1;
        break;
    case ENTER_MONITOR_RATE_MINUS:
        utils[0] = calcutil(0);
        NS_LOG_DEBUG(utils[0].util << ' ' << utils[0].avgrtt);
        monitorevent = Simulator::Schedule(MI, &ConsumerPCCMY::PCCCollect_Calc,
                                           this, ENTER_CLAC_RATE_AND_PROBE);
        monitorstate = 2;
        break;
    case ENTER_CLAC_RATE_AND_PROBE:
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
        //if ((mrate0 > 0.65)&& (gamma >0) )
        //{
        //        gamma = -(epsilon * rate);
        //        NS_LOG_DEBUG("insmarkrate larger than 0.75, gamma:" << gamma);
        //}

        // if(marksum!=0 )
        // {
        //     fincreaseflag=0;
        //     contk =1;
        // }
        // if(marksum==0 && contk>=3)
        // {
        //     gamma= std::max(gamma,rate*epsilon);
        //     NS_LOG_DEBUG("fast increase rate:");
        //     contk =1;
        // }
        // if(marksum==0 && contk >=2 && gamma>0)
        // {
        //     fincreaseflag++;
        //     NS_LOG_DEBUG("fast increase rate,fincreaseflagr:"<<fincreaseflag);
        //    // m_frequency = m_frequency + 20 * contk;
        //     gamma= std::max(gamma,20.0 * fincreaseflag);
        //     if(gamma>320)
        //     {
        //         gamma =320;
        //     }
        //     // contk =1;
        // }
        // if(mrate0 >0.9)
        // {
        //     gamma= std::min(gamma,-rate*0.2);
        //     NS_LOG_DEBUG("fast decrease rate:");
        //     contk =1;
        // }

        rate += gamma;
        NS_LOG_DEBUG("newrate: " << rate);
        if (rate < 100.0)
            rate = 100.0;
        probeevent = Simulator::ScheduleNow(&ConsumerPCCMY::PCCProbe, this,
                                            ENTER_PROBE_RATE_PLUS);
        monitorstate = 0;
        break;
    default:
        break;
    }
    m_frequency = rate;
}

void ConsumerPCCMY::ScheduleNextPacket() {
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

void ConsumerPCCMY::WillSendOutInterest(uint32_t sequenceNumber) {
    // NS_LOG_DEBUG("sendout");
    sqrMI[probestate].emplace(sequenceNumber);
    sendnum[probestate]++;
    // NS_LOG_DEBUG("emplace " << probestate - 1 << " seq " <<
    // sequenceNumber<< " done");
    ConsumerCbr::WillSendOutInterest(sequenceNumber);
}

int64_t ConsumerPCCMY::getSrcofData(shared_ptr<const Data> data) {
    auto blk =
        data->getMetaInfo().findAppMetaInfo(::ndn::tlv::AppPrivateBlock1);
    if (blk == nullptr)
        return -1;
    else
        return ::ndn::readNonNegativeInteger(*blk);
}

void ConsumerPCCMY::OnTimeout(uint32_t sequenceNum)
{
    int index = monitorstate;
 
    
    m_retxSeqs.insert(sequenceNum);
    ScheduleNextPacket();
    
}

void ConsumerPCCMY::OnData(shared_ptr<const Data> data) {
    uint64_t congmark = data->getCongestionMark();
    int index = monitorstate;
    receivenum++;
    uint32_t seqNum = data->getName().get(-1).toSequenceNumber();
    int64_t respondid = getSrcofData(data);
    
    datainfo[index].emplace(seqNum,
                            datainfo_t{seqNum, respondid, congmark, 0, 0, 0});
     ConsumerCbr::OnData(data);
   //NS_LOG_DEBUG("Received congestion mark , " << congmark);
    if (SSflag)
    {

         NS_LOG_DEBUG("enter slow start: ");
         int64_t cur_time1 = Simulator::Now().GetNanoSeconds();
         if (congmark >2 && (cur_time1 - dectime > 2.5*600000000))
          {
                //NS_LOG_DEBUG("cur_time1: " << cur_time1, << "dectime" << dectime);
                for (int j=0; j< 10; j++)
                     NS_LOG_DEBUG("SS over: ");
               SSflag = 0;
               SStime = Simulator::Now().GetNanoSeconds();
               rate = m_frequency*0.8;
               probeevent = Simulator::ScheduleNow(&ConsumerPCCMY::PCCProbe, this, 1);
          }
         else
         {
            // for (int j=0; j< 10; j++)
              //   NS_LOG_DEBUG("enter SS: ");
             
	     //
             //m_frequency=m_frequency+5;
             
             NS_LOG_DEBUG("mark: " << congmark);
             //NS_LOG_DEBUG("SSflag: "<<SSflag);
         }
    }
   
    if(!SSflag)
   {
    
    if (congmark > 3) {
           int index = monitorstate;
        //if (sqrMI[index].find(seqNum) != sqrMI[index].end()) {
             auto it = datainfo[index].find(seqNum);
            if (it != datainfo[index].end()) {
                it->second.congmark = data->getCongestionMark();
                if (data->getCongestionMark() > 30) {
                    lossnum++;
                }
                
                // NS_LOG_DEBUG("marknum: " << marknum);
                marksum += congmark;
                // }
            }
       // }
    }
    
     int64_t cur_time = Simulator::Now().GetNanoSeconds();
    
    if (datainfo[index][seqNum].rtt < 35000000)
    {
        //NS_LOG_DEBUG("datainfo[index][seqNum].rtt:" << datainfo[index][seqNum].rtt);
        varytest[0].push_back(congmark);
    }
    else if (datainfo[index][seqNum].rtt > 35000000)
    {
        varytest[1].push_back(congmark);
    }
     int num = 100;

     if (varytest[1].size() > num)
     {
	 heavycong = 0;
	 for(int i =varytest[1].size(); i>varytest[1].size()-num; i--)
	{
   	      if (varytest[1][i]>28)
	      {
		   heavycong++;
	      }
	}
     }
     
	insmarkrate = double(heavycong)/num;
      //NS_LOG_DEBUG("varytest[0]:" << varytest[0].size() << "varytest[1]: "<<varytest[1].size());
      //NS_LOG_DEBUG("calculate insmarkrate---heavyconge:" << heavycong << "insmarkrate: "<<insmarkrate);
      //NS_LOG_DEBUG("cur_time-SStime:" << cur_time-SStime);
     if (cur_time - SStime >  2000000000 && insmarkrate > 1 && cur_time - dectime > 3*600000000)
    {
         for (int j=0; j<10; j++)
            NS_LOG_DEBUG("new network state, decrease rate " );
	//NS_LOG_DEBUG("cur_time: " << cur_time << " SStime: " << SStime);
       
	       //SSflag = 1;
             //if()
	     
            // return{0,0};
            probeevent.Cancel();
            monitorevent.Cancel();
            rate = 0;
            epsilon = 0.05;
      	    theta0 = 1;
      	    omega0 = 0.05;
            delta = 0.1;
            omegaboundk = 0;
            probestate = 0;
            monitorstate = 0;
            lastdiffpositive = false;
            contk = 1;
          //  m_frequency = 100;
	    
	   rate = m_frequency*0.8;
	    m_frequency = m_frequency * 0.8;
        decflag =0;
             dectime = Simulator::Now().GetNanoSeconds();
            // Simulator::Schedule(Seconds(0), &ConsumerPCCMY::ssrate, this );
         probeevent = Simulator::Schedule(2*MI, &ConsumerPCCMY::PCCProbe, this,
                                         ENTER_PROBE_RATE_PLUS);
    }
    
    }
    
    //auto da = datainfo[index].find(seqno);

    
     
}

void ConsumerPCCMY::getdelay(Ptr<App> app, uint32_t seqno, Time delay,
                             uint32_t retxCount, int32_t hopCount) {
    int index = monitorstate;
    // NS_LOG_DEBUG("try seq " << seqno << " in " << index);
    auto it = datainfo[index].find(seqno);
    
    if (it != datainfo[index].end()) {
        it->second.rtt = delay.GetNanoSeconds();
        //NS_LOG_DEBUG("decide source---it->second.rtt_inc_rto " << it->second.rtt);
        if (it->second.rtt > 60000000) {
            it->second.src = 0;
        } else if (it->second.rtt > 30000000 && it->second.rtt < 60000000) {
            it->second.src = 1;
        }
        else if (it->second.rtt < 30000000 ) {
            it->second.src = 2;
        }
        it->second.timestamp = Simulator::Now().GetNanoSeconds();
        if (sqrMI[index].find(seqno) != sqrMI[index].end()) {
            sqrMI[index].erase(seqno);
            it->second.is_in_sqrMI = true;
        }
    }
}

void ConsumerPCCMY::SetRandomPrefix(std::string random_prefix) {
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

}  // namespace ndn
}  // namespace ns3
