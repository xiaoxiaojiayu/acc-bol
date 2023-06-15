#include "pcc5stage-merge-yjy.hpp"

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

NS_LOG_COMPONENT_DEFINE("ndn.PCC5Stage");

namespace ns3 {
namespace ndn {

::std::ostream &operator<<(::std::ostream &os, PCCProbe_MODE mode) {
    static const std::string modestr[5] = {"PL0", "PLUS", "PL1", "MINUS",
                                           "PL2"};
    os << modestr[mode];
    return os;
}

::std::ostream &operator<<(::std::ostream &os, PCCCollect_Calc_MODE mode) {
    static const std::string modestr[4] = {"NONE", "PLUS", "MINUS", "CALC"};
    os << modestr[mode];
    return os;
}

NS_OBJECT_ENSURE_REGISTERED(ConsumerPCC5Stage);

TypeId ConsumerPCC5Stage::GetTypeId() {
    static TypeId tid =
        TypeId("ns3::ndn::ConsumerPCC5Stage")
            .SetGroupName("Ndn")
            .SetParent<ConsumerCbr>()
            .AddConstructor<ConsumerPCC5Stage>()
            // .AddAttribute("rate", "", DoubleValue(200.0),
            //               MakeDoubleAccessor(&ConsumerPCC5Stage::setrate),
            //               MakeDoubleChecker<double>())
            .AddAttribute("epsilon", "", DoubleValue(0.05),
                          MakeDoubleAccessor(&ConsumerPCC5Stage::epsilon),
                          MakeDoubleChecker<double>())
            .AddAttribute(
                "RandomPrefix", "", StringValue(""),
                MakeStringAccessor(&ConsumerPCC5Stage::SetRandomPrefix),
                MakeStringChecker());

    return tid;
}

ConsumerPCC5Stage::ConsumerPCC5Stage()
    : ConsumerCbr(),
      rate(0),
      epsilon(0.05),
      theta0(1),
      omega0(0.05),
      delta(0.1),
      omegaboundk(0),
      probe_table_index(0),
      monitor_table_index(0),
      lastdiffpositive(false),
      contk(1),
      SSflag(1),
      ssfirst(false),
      SStime(0),
      dectime(0),
      insmarkrate(0),
      heavycong(0),
      decflag(1),
      receivenum(0),
      lossrate(0),
      randomloss(0),
      v(150),
      eta(0.95),
              congestionminus(0),
              congestionplus(0),
              targetqueue(15),
              // modipart(0),
              // speederflag(false),
      firsttimeset(2),
      mrate0(0),
              beta(0.08),
              srcnum(3)
        {
      
    m_firstInterestDataDelay.ConnectWithoutContext(
        MakeCallback(&ConsumerPCC5Stage::getdelay_inc_rto, this));
    m_lastRetransmittedInterestDataDelay.ConnectWithoutContext(
        MakeCallback(&ConsumerPCC5Stage::getdelay_no_rto, this));
    maxrtt = MilliSeconds(50);
    minrtt = MilliSeconds(50);
}

ConsumerPCC5Stage::~ConsumerPCC5Stage() {}

void ConsumerPCC5Stage::StartApplication() {
    Simulator::ScheduleNow(&ConsumerPCC5Stage::ssrate, this);
    // probeevent = Simulator::ScheduleNow(&ConsumerPCC5Stage::PCCProbe, this,
    //                                     ENTER_PROBE_RATE_PAYLOAD_0);
    ConsumerCbr::StartApplication();
}

void ConsumerPCC5Stage::resettable() {
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
    lossrate = 0;
    mrate0 = 0;
    queueold_inc = new int[srcnum];
    queueold_rec= new int[srcnum];
    queuestart_inc=new int[srcnum];
    queuestart_rec=new int[srcnum];
}

void ConsumerPCC5Stage::updateMaxMinRttFromDatainfo(int index) {
    if (datainfo[index].empty()) {
        NS_LOG_DEBUG("datainfo[info] is empty, can't update maxminRtt");
        return;
    }
    int64_t curmin = (~0UL >> 1);
    int64_t curmax = 0;
    for (auto i : datainfo[index]) {
        if (i.second.rtt < curmin) {
            curmin = i.second.rtt;
        }
        if (i.second.rtt > curmax) {
            curmax = i.second.rtt;
        }
    }
    //maxrtt = NanoSeconds(70000000);
    //minrtt = NanoSeconds(20000000);
    NS_LOG_DEBUG("maxminRtt updated, max: " << maxrtt << ", min: " << minrtt);
}

void ConsumerPCC5Stage::updateMI() {
    MI = NanoSeconds(maxrtt.GetNanoSeconds());
    rttdiff = NanoSeconds((maxrtt - minrtt).GetNanoSeconds());
    NS_LOG_DEBUG("MI updated, MI: " << MI << " rttdiff: " << rttdiff);
}

void ConsumerPCC5Stage::PCCProbe(PCCProbe_MODE mode) {
    NS_LOG_DEBUG("Probe mode: " << mode);
    switch (mode) {
    case ENTER_PROBE_RATE_PAYLOAD_0:
        m_frequency = rate;
        probeevent = Simulator::Schedule(rttdiff, &ConsumerPCC5Stage::PCCProbe,
                                         this, ENTER_PROBE_RATE_PLUS);
        probe_table_index = 0;
        break;
    case ENTER_PROBE_RATE_PLUS:
        m_frequency = rate * (1 + epsilon);
        probeevent = Simulator::Schedule(MI, &ConsumerPCC5Stage::PCCProbe, this,
                                         ENTER_PROBE_RATE_PAYLOAD_1);
        probe_table_index = 1;

        resettable();
        monitorevent =
            Simulator::Schedule(minrtt, &ConsumerPCC5Stage::PCCCollect_Calc,
                                this, ENTER_MONITOR_RATE_PLUS);
        break;
    case ENTER_PROBE_RATE_PAYLOAD_1:
        m_frequency = rate;  // * (1 + epsilon);//rate;
        probeevent = Simulator::Schedule(rttdiff, &ConsumerPCC5Stage::PCCProbe,
                                         this, ENTER_PROBE_RATE_MINUS);
        probe_table_index = 0;

        break;
    case ENTER_PROBE_RATE_MINUS:
        m_frequency = rate * (1 - epsilon);
        probeevent = Simulator::Schedule(MI, &ConsumerPCC5Stage::PCCProbe, this,
                                         ENTER_PROBE_RATE_PAYLOAD_2);
        monitorevent =
            Simulator::Schedule(minrtt, &ConsumerPCC5Stage::PCCCollect_Calc,
                                this, ENTER_MONITOR_RATE_MINUS);

        probe_table_index = 2;
        break;
    case ENTER_PROBE_RATE_PAYLOAD_2:
        m_frequency = rate;  //* (1 - epsilon);//rate;
        probe_table_index = 0;

        break;
    //default:
        break;
    }
    NS_LOG_DEBUG("freq: " << m_frequency);
}

void ConsumerPCC5Stage::PCCCollect_Calc(PCCCollect_Calc_MODE mode) {
    NS_LOG_DEBUG("Collect_calculate mode: " << mode);

    utilret utils[2];
    double gamma = 0;
    double utildiff = 0;
    double rttdouble = 0;
    double omega = 0;
    double marksum_old = 0;
            double modipart = 0;

    switch (mode) {
    case ENTER_MONITOR_NONE:
        monitor_table_index = 0;
        break;
    case ENTER_MONITOR_RATE_PLUS:
        monitorevent = Simulator::Schedule(MI + rttdiff,
                                           &ConsumerPCC5Stage::PCCCollect_Calc,
                                           this, ENTER_MONITOR_NONE);
        monitor_table_index = 1;
        break;
    case ENTER_MONITOR_RATE_MINUS:
        updateMaxMinRttFromDatainfo(1);
        updateMI();
        monitorevent = Simulator::Schedule(MI + rttdiff,
                                           &ConsumerPCC5Stage::PCCCollect_Calc,
                                           this, ENTER_CLAC_RATE_AND_PROBE);
        monitor_table_index = 2;
        break;
    case ENTER_CLAC_RATE_AND_PROBE:
        updateMaxMinRttFromDatainfo(2);
        updateMI();
        utils[0] = calcutil(0);
        utils[1] = calcutil(1);
        NS_LOG_DEBUG("util[r+]=" << utils[0].util
                                 << " util[r-]=" << utils[1].util);

        // if (utils[0].lossrate > 0.1 && utils[1].lossrate > 0.1) {
        //     NS_LOG_DEBUG("lossrate[+][-] both exceed 0.1, discard all states");
        //     rate = 0.5 * rate;
        //     if (rate < 100.0)
        //         rate = 100.0;
        //     NS_LOG_DEBUG("new rate: " << rate);
        //     contk = 1;
        //     probeevent = Simulator::ScheduleNow(
        //         &ConsumerPCC5Stage::PCCProbe, this, ENTER_PROBE_RATE_PAYLOAD_0);
        //     monitorevent = Simulator::ScheduleNow(
        //         &ConsumerPCC5Stage::PCCCollect_Calc, this, ENTER_MONITOR_NONE);
        //     break;
        // }

                if (congestionplus > congestionminus)
                {
                    if(congestionplus < targetqueue)
                    {
                            modipart = 2 * beta * rate * (targetqueue - (congestionplus + congestionminus) / 2.0);
                    }
                    else
                    {
                          modipart = (congestionplus-targetqueue)* beta * rate * (targetqueue - (congestionplus + congestionminus) / 2.0);
                    }
                    
                    NS_LOG_DEBUG("Normal queue, add modifypart, congestionplus: " << congestionplus << "congestionminus: " << congestionminus << " modipart: " << modipart);
                }
                else
                {
                    if (congestionminus < targetqueue)
                    {
                        modipart = 2 * beta * rate * (congestionplus - congestionminus);
                    }
                    else
                    {
                        modipart = 2 * beta * rate * (congestionplus - congestionminus) + (congestionminus-targetqueue) * beta * rate * (targetqueue - (congestionminus + congestionplus) / 2.0);
                    }

                    NS_LOG_DEBUG("Normal queue, add modifypart, congestionplus: " << congestionplus << "congestionminus: " << congestionminus << " modipart: " << modipart);
                }

                utildiff = utils[0].util - utils[1].util + modipart;

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

        // if (mrate0 > 0.98) {
        //     gamma = -0.2 * rate;
        //     contk = 1;
        //     omegaboundk = 0;
        // }

        NS_LOG_DEBUG("gamma:" << gamma);
        NS_LOG_DEBUG("omegaboundk:" << omegaboundk);

        rate += gamma;

        if (rate < 100.0)
            rate = 100.0;
        NS_LOG_DEBUG("new rate: " << rate);

         int64_t cur_time = Simulator::Now().GetNanoSeconds();
        NS_LOG_DEBUG("curtime:"<< cur_time);
        if(cur_time>10000000000 && cur_time<13000000000 )
        {

                    if (contk >= 3 && marksum <= 20 && gamma > 0)
                    {
            SSflag = 1;
            contk = 1;
            Simulator::ScheduleNow(&ConsumerPCC5Stage::ssrate, this);
        } else {
            probeevent = Simulator::ScheduleNow(
                &ConsumerPCC5Stage::PCCProbe, this, ENTER_PROBE_RATE_PAYLOAD_0);
            monitorevent = Simulator::ScheduleNow(
                &ConsumerPCC5Stage::PCCCollect_Calc, this, ENTER_MONITOR_NONE);
        }
        }
        else{
            probeevent = Simulator::ScheduleNow(&ConsumerPCC5Stage::PCCProbe, this,
                                            ENTER_PROBE_RATE_PAYLOAD_0);
        monitorevent = Simulator::ScheduleNow(
            &ConsumerPCC5Stage::PCCCollect_Calc, this, ENTER_MONITOR_NONE);
        }
        break;
    //default:
        break;
    }
    m_frequency = rate;
}

double ConsumerPCC5Stage::funcA(double rate, double drtt_dt1, double weight1,
                                double drtt_dt2, double weight2,
                                double markrate) {
    NS_LOG_DEBUG(rate << " " << drtt_dt1 << " " << drtt_dt2 << " " << markrate);
    double w = 0;
    double maxcong = 0;
    if (drtt_dt1 > drtt_dt2) {
        maxcong = drtt_dt1;
        w = weight1;
        return pow(rate, 0.9) - 900 * rate * weight1 * drtt_dt1 -
               11.35 * markrate * rate * weight1;
    } else {
        maxcong = drtt_dt2;
        w = weight2;
        return pow(rate, 0.9) - 900 * rate * weight2 * drtt_dt2 -
               11.35 * markrate * rate * weight2;
    }

    // return pow(rate, 0.9) - 900 * rate * weight2 * drtt_dt2 -
    //        11.35 * markrate * rate * weight2;
}

double ConsumerPCC5Stage::func1(double rate, double markrate, double lossrate) {
    // NS_LOG_DEBUG(rate << " " << drtt_dt1 << " " << drtt_dt2 << " " <<
    // markrate);
    return pow(rate, 0.9) - 900 * markrate * rate - 11.35 * lossrate * rate;
}

double ConsumerPCC5Stage::funcB(double rate, double drtt_dt1, double weight1,
                                double drtt_dt2, double weight2,
                                double lossrate1,double lossrate2) {
    NS_LOG_DEBUG("rate : "<<  rate << " avgmark1 : " << drtt_dt1 << " avgmark2 :" << drtt_dt2 << "mrate[0]:  " << lossrate1 << "mrate[1]:  " << lossrate2);
    // return 1.2 * rate * log(rate) - 134 * rate * (weight1 * drtt_dt1 +
    // weight2 * drtt_dt2) - 121.8 * lossrate * rate;  // ;//
    // return log(rate)/13.0 - 900 * rate * (weight1 * drtt_dt1 + weight2 *
    // drtt_dt2) - 11.35 * lossrate * rate;
    //if (mrate0 > 0.5)
    //    return -113.5 * lossrate * rate;
    //else
         return 15 * pow(rate, 0.9) -
            30 * rate * (weight1 * drtt_dt1 + weight2 * drtt_dt2) -30 * rate*lossrate1;//(lossrate1 * weight1+ lossrate2*weight2);
}  // 0.1 * rate*log(rate)

double ConsumerPCC5Stage::funcC(double rate, double drtt_dt1, double weight1,
                                double drtt_dt2, double weight2, double lrate) {
    NS_LOG_DEBUG(rate << " " << drtt_dt1 << " " << drtt_dt2 << " " << lrate);
    return 2.5 * rate - 100 * rate * (weight1 * drtt_dt1 + weight2 * drtt_dt2) -
           11.35 * lrate * rate;
}

void ConsumerPCC5Stage::ssrate() {
    if (SSflag) {
        if(ssfirst)
        {
        NS_LOG_DEBUG("enter SSrate: ");
        int deltarate = 10 * (1 << SSflag / 2);
        if (deltarate > 320)
            deltarate = 320;
        else if (deltarate < 0)
            deltarate = 0;
        NS_LOG_DEBUG("m_frequency:" << m_frequency << " + " << deltarate);
        //m_frequency = m_frequency + deltarate;
         m_frequency = m_frequency*(1+epsilon);
        }
        if(!ssfirst)
        {
            m_frequency = m_frequency*(1+epsilon);
        }
        
        NS_LOG_DEBUG("m_frequency: " << m_frequency);
        SSflag++;
        ScheduleNextPacket();
        Simulator::Schedule(MI, &ConsumerPCC5Stage::ssrate, this);
        updateMaxMinRttFromDatainfo(0);
        updateMI();
    }
}

double ConsumerPCC5Stage::m(const double &diff) {
    // if (lastdiffpositive ^ (diff > 0)) {
    //     contk = 1;
    //     omegaboundk = 0;
    // } else {
    //     contk++;
    // }

    // lastdiffpositive = (diff > 0);

    // NS_LOG_DEBUG("contk:" << contk);

    // if (diff < 0) {
    //     if (contk <= 3)
    //         return diff;
    //     else
    //         return (0.5 * contk - 1.5) * diff;
    // } else {
    //     if (contk <= 3)
    //         return contk * diff;
    //     else
    //         return (2.5 * contk - 3) * diff;
    // }
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

utilret ConsumerPCC5Stage::calcutil(int i) {
    NS_LOG_DEBUG("entercalcutil: ");
    i++;

    const size_t rttnum = datainfo[i].size();
    double predictloss =0;
    lossrate = (double)sqrMI[i].size() / (double)(sendnum[i] + 1);
    // double(lossnum) / sendnum[i];  //
    if (rttnum < 10 || lossrate < 0.02)
        lossrate = 0;
    const double markrate = double(marknum) / receivenum;  //(sendnum[i]+1);

    NS_LOG_DEBUG("rttnum: " << rttnum);
    NS_LOG_DEBUG("sqrMI[i].size(): " << (double)sqrMI[i].size()
                                     << " sendnum[i]: " << sendnum[i]);
    NS_LOG_DEBUG("lossrate: " << lossrate);
    NS_LOG_DEBUG("marknum: " << marknum << " receivenum: " << receivenum);
    NS_LOG_DEBUG("markrate: " << markrate);

    double nullval;
    int rttnum1 = 0;
    int rttnum2 = 0;
    int64_t rttsum1 = 0;
    int64_t rttsum2 = 0;
    int64_t rttsum = 0;
    double weight1 = 0;
    double weight2 = 0;
    double rttgrad1 = 0;
    double rttgrad2 = 0;
    double avgmark1 = 0, avgmark2 = 0;
    double heavymark1=0, heavymark2=0;
    double lossestimate = 0;
            int groupsize = 10;
            int tmpindex = 0;

    //const size_t srcnum = 3;
    
    double *avgrtt = new double[srcnum];
    int *recvnum = new int[srcnum];
    int *countnum = new int[srcnum];
    double *weight = new double[srcnum];
    int *mknum = new int[srcnum];
    int *mksum = new int[srcnum];
    double *mrate = new double[srcnum];
    double *smrtt = new double[srcnum];
    double *maxqueue = new double[srcnum];
            double *totalmark = new double[srcnum];
            double *avgcontqueue = new double[100]{0};
            int *countindex = new int[srcnum];


    for (int src = 0; src < srcnum; src++) {
        avgrtt[src] = 0;
        weight[src] = 0;
        mrate[src] = 0;
        recvnum[src] = 0;
        mksum[src] = 0;
        mknum[src] = 0;
        smrtt[src] = 0;
                maxqueue[src] = 0;
                countnum[src] = 0;
                totalmark[src] = 0;

                avgcontqueue[src] = 0;
                countindex[src] = 0;
    }

    double *x1 = new double[rttnum];
    double *y1 = new double[rttnum];
    double *x2 = new double[rttnum];
    double *y2 = new double[rttnum];

    auto it = datainfo[i].cbegin();
    if(rttnum > 5)
    {
        maxrtt = MilliSeconds(5);
        minrtt = MilliSeconds(90);
    }
    
    while (it != datainfo[i].cend()) {
        if(NanoSeconds(it->second.rtt) > maxrtt && it->second.rtt < 200000000)
        { maxrtt = NanoSeconds(it->second.rtt);
        }
        if(NanoSeconds(it->second.rtt) < minrtt)        {
            minrtt = NanoSeconds(it->second.rtt);
        }
        int cursrc = it->second.src;
                countindex[cursrc] += 1;
                avgcontqueue[cursrc] += it->second.congmark;
                if (countindex[cursrc] > groupsize)
                {

                    tmpindex = cursrc + groupsize;
                    if (avgcontqueue[cursrc] / groupsize > avgcontqueue[tmpindex])
                    {
                        avgcontqueue[tmpindex] = avgcontqueue[cursrc] / double(groupsize);
                        NS_LOG_DEBUG("findnew avgcontqueue: srcnum: " << cursrc << "avgcontqueue[tmpindex]: " << avgcontqueue[tmpindex]);
                    }

                    avgcontqueue[cursrc] = 0;
                    countindex[cursrc] = 0;
                }
        smrtt[cursrc] += it->second.rtt / 1000000.0;
       
         if (it->second.congmark >30) {
            
            lossestimate+= 1;
         }
        if (it->second.congmark >= 10) {
            mknum[cursrc] += 1;
            mksum[cursrc] += it->second.congmark;
            
        }
        recvnum[cursrc] += 1;
        countnum[cursrc]+=1;
        
        // NS_LOG_DEBUG("enter while circle " );
        int64_t rtt = it->second.rtt;
        int64_t ts = it->second.timestamp;
        if (cursrc == 0) {
            // NS_LOG_DEBUG("The source is r0: " );
            x1[rttnum1] = ts / 1000000000.0;
            // y1[rttnum1] = rtt / 1000000000.0;
            y1[rttnum1] = it->second.congmark;
            if (it->second.congmark >30) {
            
            heavymark1+= 1;
            
        }
        if(it->second.congmark>maxqueue[cursrc])
                maxqueue[cursrc]=it->second.congmark;
                
            rttsum1 += rtt;
            rttnum1++;
            rttsum += rtt;
        } else if (cursrc == 1) {
            // NS_LOG_DEBUG("The source is p0: " );

            x2[rttnum2] = ts / 1000000000.0;
            // y2[rttnum2] = rtt / 1000000000.0;
            y2[rttnum2] = it->second.congmark;
            if (it->second.congmark >23) {
            
            heavymark2+= 1;
        }
        if(it->second.congmark>maxqueue[cursrc])
                maxqueue[cursrc]=it->second.congmark;

            rttsum2 += rtt;
            rttnum2++;
            rttsum += rtt;
        }else if (cursrc == 2) {
                if(it->second.congmark>maxqueue[cursrc])
                    maxqueue[cursrc]=it->second.congmark;

                                                
        }

        it++;
    }
            double tmp = 0;
            for (int src = 0; src < srcnum; src++)
            {
        avgrtt[src] = smrtt[src] / (1 + recvnum[src]);
        weight[src] = double(countnum[src]) / (1 + rttnum);
        mrate[src] = double(mksum[src]) / (1 + recvnum[src]);
                mrate[src] = mrate[src] / 32.0;
                tmpindex = groupsize + src;
                tmp += weight[src] * avgcontqueue[tmpindex];
                NS_LOG_DEBUG("srcnum: " << src << "avgcontqueue[tmpindex]: " << avgcontqueue[tmpindex]);
                avgcontqueue[src] = avgcontqueue[tmpindex];
            }

            for (int src = 0; src < srcnum; src++)
            {

                // avgrtt[src] = smrtt[src] / (1 + recvnum[src]);
                // weight[src] = double(countnum[src]) / (1 + rttnum);
                // mrate[src] = double(mksum[src]) / (1 + recvnum[src]);
                // mrate[src]=mrate[src]/15.0;
            }
            if (i == 1)
            {
                congestionplus = tmp;
                NS_LOG_DEBUG("congestionplus: " << congestionplus);
            }
            else
            {
                congestionminus = tmp;
                NS_LOG_DEBUG("congestionminus: " << congestionminus);
    }
    //lossrate = lossestimate/(recvnum[1]+recvnum[2]+recvnum[3]);
    std::stringstream ssx1, ssx2;
    std::stringstream ssy1, ssy2;
    for (int j = 0; j < rttnum1; j++) {
        ssx1 << x1[j] << '\t';
        ssy1 << y1[j] << '\t';
    }
    for (int j = 0; j < rttnum2; j++) {
        ssx2 << x2[j] << '\t';
        ssy2 << y2[j] << '\t';
    }

    NS_LOG_DEBUG("ssx1.str()"<<ssx1.str());
    NS_LOG_DEBUG("ssy1.str()"<<ssy1.str());
    NS_LOG_DEBUG("ssx2.str()"<<ssx2.str());
    NS_LOG_DEBUG("ss.sy2tr()"<<ssy2.str());
    
        for(int j=0; j<rttnum1; j++)
        {
            y1[0]+=y1[j];
            //NS_LOG_DEBUG("y1[j]: " << y1[j]);
        }
        avgmark1 = double(y1[0])/(32.0*(rttnum1+1));
        NS_LOG_DEBUG("avgmark1: " << avgmark1<<"y1[0]:"<<y1[0]<<"rttnum1:"<<rttnum1);

         for(int j=0; j<rttnum1; j++)
        {
            y2[0]+=y2[j];
            //NS_LOG_DEBUG("y1[j]: " << y2[j]);
        }

        avgmark2 = double(y2[0])/(32.0*(rttnum2+1));
        NS_LOG_DEBUG("avgmark2: " << avgmark2<<"y2[0]:"<<y2[0]<<"rttnum2:"<<rttnum2);
        

    gsl_fit_linear(x1, 1, y1, 1, rttnum1, &nullval, &rttgrad1, &nullval,
                   &nullval, &nullval, &nullval);
    gsl_fit_linear(x2, 1, y2, 1, rttnum2, &nullval, &rttgrad2, &nullval,
                   &nullval, &nullval, &nullval);

    delete[] x1;
    delete[] y1;
    delete[] x2;
    delete[] y2;

    if (abs(rttgrad1) < 100 || isnan(rttgrad1))
        rttgrad1 = 0;
    if (abs(rttgrad2) < 100 || isnan(rttgrad2))
        rttgrad2 = 0;

    weight1 = double(rttnum1) / (rttnum + 1);
    weight2 = double(rttnum2) / (rttnum + 1);

    NS_LOG_DEBUG("rttsum1: " << rttsum1 << " rttnum1: " << rttnum1
                             << " weight1: " << weight1
                                     << " markgrad1: " << rttgrad1 << "maxqueue[0]: " << maxqueue[0]);
    NS_LOG_DEBUG("rttsum2: " << rttsum2 << " rttnum2: " << rttnum2
                             << " weight2: " << weight2
                                     << " markgrad2: " << rttgrad2 << "maxqueue[1]: " << maxqueue[1]);
    NS_LOG_DEBUG("avgmark1: " << avgmark1 << " recenum[0]: " << recvnum[0]
                              << " weight1: " << weight[0] << " mknum[0]: "
                              << mknum[0] << " mkrate[0]: " << mrate[0]);
    NS_LOG_DEBUG("avgmark2: " << avgmark2 << " recenum[1]: " << recvnum[1]
                              << " weight2: " << weight[1] << " mknum: "
                              << mknum[1] << " mkrate[1]: " << mrate[1]);

    i--;
    mrate0 = std::max(heavymark1/(double)rttnum1, heavymark2/(double)rttnum2);
    NS_LOG_DEBUG("mrate0: " << mrate0);
    /*if(mrate0 < 0.001)
    {
        if(randomloss!=0 && randomloss>lossrate)
            randomloss=lossrate;
        if(randomloss==0)
            randomloss=lossrate;

    }*/

    // NS_LOG_DEBUG("randomloss: " << randomloss);
    //  return {funcA(rate * (1 + (2 * i - 1) * epsilon), rttgrad1, weight1,
    //               rttgrad2, weight2, markrate),
    //          (double)rttsum / rttnum};
    
    return {funcB(rate * (1 + (1 - 2 * i) * epsilon), mrate[0], weight[1],
                  mrate[1], weight[2], lossrate,mrate[1])};//mrate[0]
    // return {funcC(rate * (1 + (2 * i - 1) * epsilon), rttgrad1, weight1,
    //              rttgrad2, weight2, lossrate)};
    //   (double)rttsum / rttnum};

    // return {func1(rate * (1 + (2 * i - 1) * epsilon), markrate, lossrate),
    //         (double)rttsum / rttnum};


}
void ConsumerPCC5Stage::ScheduleNextPacket() {
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

void ConsumerPCC5Stage::WillSendOutInterest(uint32_t sequenceNumber) {
    // NS_LOG_DEBUG("sendout");
    sqrMI[probe_table_index].emplace(sequenceNumber);
    sendnum[probe_table_index]++;
    // NS_LOG_DEBUG("emplace " << probestate - 1 << " seq " <<
    // sequenceNumber<< " done");
    ConsumerCbr::WillSendOutInterest(sequenceNumber);
}

int64_t ConsumerPCC5Stage::getSrcofData(shared_ptr<const Data> data) {
    auto blk =
        data->getMetaInfo().findAppMetaInfo(::ndn::tlv::AppPrivateBlock1);
    if (blk == nullptr)
        return -1;
    else
        return ::ndn::readNonNegativeInteger(*blk);
}

void ConsumerPCC5Stage::OnData(shared_ptr<const Data> data) {
    int index = monitor_table_index;
    receivenum++;

    uint64_t congmark = data->getCongestionMark();
    uint32_t seqNum = data->getName().get(-1).toSequenceNumber();
    int64_t respondid = getSrcofData(data);
    //NS_LOG_DEBUG("respondid---: "<<respondid);
    datainfo[index].emplace(seqNum,
                            datainfo_t{seqNum, respondid, congmark, 0, 0, 0});
    ConsumerCbr::OnData(data);

    if (SSflag) {
        // NS_LOG_DEBUG("In SS: ");
        // NS_LOG_DEBUG("Received congestion mark , " << congmark);
        // int64_t cur_time1 = Simulator::Now().GetNanoSeconds();
        if (!ssfirst) {  //&& (cur_time1 - dectime > 2.5 * 600000000)
            // NS_LOG_DEBUG("cur_time1: " << cur_time1, << "dectime" <<
            // dectime);
                    /// if (!ssfirst || SSflag > 20) {
                    if (congmark > 7)
                {
                for (int j = 0; j < 10; j++)
                    NS_LOG_DEBUG("SS over: ");

                // SStime = Simulator::Now().GetNanoSeconds();
                int deltarate = 10 * (1 << SSflag / 2);
                if (deltarate > 320)
                    deltarate = 320;
                else if (deltarate < 0)
                    deltarate = 0;
                NS_LOG_DEBUG("withdraw last delta: " << deltarate);
                rate = m_frequency * 0.9;
                SSflag = 0;
                probeevent =
                    Simulator::ScheduleNow(&ConsumerPCC5Stage::PCCProbe, this,
                                           ENTER_PROBE_RATE_PAYLOAD_0);
               
            }
        }
         if(ssfirst)
         {
            if(  congmark > 3)
            {
               
            for (int j = 0; j < 10; j++)
                    NS_LOG_DEBUG("SS over: ");

                // SStime = Simulator::Now().GetNanoSeconds();
                int deltarate = 10 * (1 << SSflag / 2);
                if (deltarate > 320)
                    deltarate = 320;
                else if (deltarate < 0)
                    deltarate = 0;
                NS_LOG_DEBUG("withdraw last delta: " << deltarate);
                rate = m_frequency * 0.9;
                SSflag = 0;
                probeevent =
                    Simulator::ScheduleNow(&ConsumerPCC5Stage::PCCProbe, this,
                                           ENTER_PROBE_RATE_PAYLOAD_0);
                ssfirst = false;
            }
         }
         
    } else {
        // marksum += congmark;
        int index = monitor_table_index;
        // if(queuestart_inc==1 && index==1)
        // {
        //     queueold_inc = congmark;
        //     queuestart_inc = 0;
        // }
        // if(queuestart_rec==1 && index==2)
        // {
        //     queueold_rec = congmark;
        //     queuestart_rec = 0;
        // }
        if (congmark > 1) {
            int index = monitor_table_index;
            // NS_LOG_DEBUG("try seq " << seqno << " in " << index);
            auto it = datainfo[index].find(seqNum);
            if (it != datainfo[index].end()) {
                it->second.congmark = data->getCongestionMark();
                if (data->getCongestionMark() > 30) {
                    lossnum++;
                }
                marknum++;
                // NS_LOG_DEBUG("marknum: " << marknum);
                marksum += congmark;
                // }
            }
        }

        int64_t cur_time = Simulator::Now().GetNanoSeconds();

        if (datainfo[index][seqNum].rtt < 35000000) {
            // NS_LOG_DEBUG("datainfo[index][seqNum].rtt:" <<
            // datainfo[index][seqNum].rtt);
            varytest[0].push_back(congmark);
        } else if (datainfo[index][seqNum].rtt > 35000000) {
            varytest[1].push_back(congmark);
        }
        int num = 50;

        if (varytest[1].size() > num) {
            heavycong = 0;
            for (int i = varytest[1].size(); i > varytest[1].size() - num;
                 i--) {
                if (varytest[1][i] > 30) {
                    heavycong++;
                }
            }
        }

        insmarkrate = double(heavycong) / num;
        //NS_LOG_DEBUG("varytest[0]:" << varytest[0].size() << "varytest[1]:"<<varytest[1].size()); 
        //NS_LOG_DEBUG("calculate insmarkrate---heavyconge:" << heavycong << "insmarkrate:"<<insmarkrate); 
        //NS_LOG_DEBUG("cur_time-SStime:" <<  cur_time-SStime); 
        if (cur_time - SStime > 2000000000 && insmarkrate > 1 && cur_time - dectime > 3 * 600000000) {
            for (int j = 0; j < 10; j++)
                NS_LOG_DEBUG("new network state, decrease rate ");
            // NS_LOG_DEBUG("cur_time: " << cur_time << " SStime: " <<   SStime);

            // SSflag = 1;
            // if()

            // return{0,0};
            probeevent.Cancel();
            monitorevent.Cancel();
            rate = 0;
            epsilon = 0.05;
            theta0 = 1;
            omega0 = 0.05;
            delta = 0.1;
            omegaboundk = 0;
            probe_table_index = 0;
            monitor_table_index = 0;
            lastdiffpositive = false;
            contk = 1;
            //  m_frequency = 100;

            rate = m_frequency * 0.9;
            m_frequency = m_frequency * 0.9;
            decflag = 0;
            dectime = Simulator::Now().GetNanoSeconds();
            // Simulator::Schedule(Seconds(0), &ConsumerPCCMY::ssrate, this
            //);
             probeevent =
                Simulator::Schedule(2 * MI, &ConsumerPCC5Stage::PCCProbe,
                this,
                                    ENTER_PROBE_RATE_PLUS);
        }
    }

    // auto da = datainfo[index].find(seqno);
}

void ConsumerPCC5Stage::getdelay_inc_rto(Ptr<App> app, uint32_t seqno,
                                         Time delay, uint32_t retxCount,
                                         int32_t hopCount) {
    int index = monitor_table_index;
    // NS_LOG_DEBUG("try seq " << seqno << " in " << index);
    auto it = datainfo[index].find(seqno);
    it->second.rtt_inc_rto = delay.GetNanoSeconds();
    if(firsttimeset>0)
    {
        maxrtt = MilliSeconds(5);
    minrtt = MilliSeconds(90);
    if(NanoSeconds(it->second.rtt_inc_rto) > maxrtt && it->second.rtt_inc_rto < 200000000)
        {
            
            maxrtt = NanoSeconds(it->second.rtt_inc_rto);
            firsttimeset--;
            
        }   
        if(NanoSeconds(it->second.rtt_inc_rto) < minrtt)
        {
            minrtt = NanoSeconds(it->second.rtt_inc_rto);
            firsttimeset--;
        }
    }
    

    //it->second.rtt_inc_rto = delay.GetNanoSeconds();
        

    //NS_LOG_DEBUG("total data--it->second.rtt_inc_rto " << it->second.rtt_inc_rto);
    if (it != datainfo[index].end()) {
        

        //("decide source---it->second.rtt_inc_rto " << it->second.rtt_inc_rto);
        if (it->second.rtt_inc_rto > 35000000) {
            it->second.src = 0;
        } else if (it->second.rtt_inc_rto > 30000000 && it->second.rtt_inc_rto < 80000000) {
            it->second.src = 2;
        }
        else if (it->second.rtt_inc_rto < 35000000 ) {
            it->second.src = 1;
        }
        it->second.timestamp = Simulator::Now().GetNanoSeconds();
        if (sqrMI[index].find(seqno) != sqrMI[index].end()) {
            sqrMI[index].erase(seqno);
            it->second.is_in_sqrMI = true;
        }
    }
}

void ConsumerPCC5Stage::getdelay_no_rto(Ptr<App> app, uint32_t seqno,
                                        Time delay, int32_t hopCount) {
    int index = monitor_table_index;
    // NS_LOG_DEBUG("try seq " << seqno << " in " << index);
    auto it = datainfo[index].find(seqno);
    if (it != datainfo[index].end()) {
        it->second.rtt = delay.GetNanoSeconds();
    }
}

void ConsumerPCC5Stage::SetRandomPrefix(std::string random_prefix) {
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
