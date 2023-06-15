#ifndef RL_H
#define RL_H
#include "ns3/ndnSIM/model/ndn-common.hpp"
#include "ns3/ndnSIM/apps/ndn-consumer.hpp"
#include "ns3/ndnSIM/apps/ndn-consumer-window.hpp"
#include "ns3/ns3-ai-module.h"
#include "ns3/traced-value.h"
#include "ns3/watchdog.h"
#include <memory>
#include <fstream>

namespace ns3 {
    namespace ndn {
        class ConsumerCCs;

        typedef struct {
            double   cWndSum;
            double   avgDelay;
            uint32_t DataNum;
            uint32_t InflightNum;
            uint32_t NackNum;
            uint32_t TimeoutNum;
            uint32_t congesLevelSum;
            uint32_t dataSizeSum;
        }Packed DDPGParam;

        typedef struct {
            //bool doAction;
            double new_cWnd;
        }Packed DDPGAct;

        enum CCType {
            AIMD, RL,
            ECP, none
        };

        class TransParam2Py : Ns3AIRL<DDPGParam, DDPGAct> {
        public:
            TransParam2Py(uint16_t id);
            void GetActFromRLModule(DDPGAct*);
            void SendParam2RLModule(DDPGParam*);
            friend ConsumerCCs;
        };

        class ConsumerCCs : public ConsumerWindow {
        public:
            static TypeId GetTypeId();
            ConsumerCCs();
            void printCollectInfo(int no, bool dg);

        protected:
            virtual void OnData(shared_ptr<const Data> data) override;
            virtual void OnTimeout(uint32_t sequenceNum) override;
            virtual void OnNack(shared_ptr<const lp::Nack> nack) override;
            virtual void ScheduleNextPacket() override;
            virtual void WillSendOutInterest(uint32_t sequenceNumber)override;


            void WindowIncrease(std::string dataname);
            void WindowDecrease(bool nacktrig);
            void changeprefix();
            uint8_t log_mask;

        private:
            CCType m_ccAlgorithm;
            double m_beta;
            double m_ssthresh;
            //double m_addRttSuppress;
            //bool m_reactToCongestionMarks;
            //bool m_useCwa;
            uint32_t m_highData;
            double m_recPoint;

        private:
            bool adjust;
            std::set<uint32_t> nackDeSeq;

            double localFullDelay;
            double localLastDelay;
            static void localLastDelayRecordCb(Ptr<App> /* app */, uint32_t /* seqno */, Time /* delay */,
                                               int32_t /*hop count*/);
            static void localFullDelayRecordCb(Ptr<App> /* app */, uint32_t /* seqno */, Time /* delay */,
                                               uint32_t /*retx count*/, int32_t /*hop count*/);

            void SetWatchDog(double tInterval);
            void SetFrequency(double f);
            void SetRandomize(const std::string& value);
            void SetTrans(uint16_t id);
            void SetLog2File(std::string path);
            void SetRandomPrefix(std::string random_prefix);
            bool randprefix;
            std::map<double, std::string> intervalmap;
            std::vector<double> intervalp;

            std::set<uint32_t> ustc2seq;

            std::ofstream logstream;

            double frequency;
            Ptr<RandomVariableStream> randomSend;

            Watchdog cwndChangeWD;
            double watchdogt;

            DDPGParam collectInfo;
            DDPGAct action;
            std::unique_ptr<TransParam2Py> transclass;
            double alpha;
	        double flag=1;
            uint32_t n=1;
            //void dynamic_change();
            friend void cwndChangeWDCallback(ConsumerCCs*);
        };
    }
}
#endif
