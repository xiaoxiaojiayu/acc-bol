#include <queue>
#include <unordered_set>

#include "ns3/ndnSIM/apps/ndn-consumer-cbr.hpp"

namespace ns3 {
namespace ndn {

typedef double (*utilfunc)(double rate, double drtt_dt, double lrate);

class ConsumerPCC : public ConsumerCbr {
  public:
    ConsumerPCC();
    virtual ~ConsumerPCC();

    static TypeId GetTypeId();

    void setrate(double ir);

  protected:

    void changeprefix();

    // 0: r 1: r-er, 2: r+er
    void PCCProbe(uint8_t mode);

    // 0:prepare 1: r-er 2: r+er
    void PCCCollect_Calc(uint8_t mode);

    void resettable();

    void updateMI();

    void getdelay(Ptr<App> app, uint32_t seqno, Time delay, uint32_t retxCount,
                  int32_t hopCount);

    void WillSendOutInterest(uint32_t sequenceNumber) override;

    double func(double rate, double drtt_dt, double lrate);

    double m(const double &diff);

    void SetRandomPrefix(std::string random_prefix);
    bool randprefix;
    //void SetRandomize(const std::string& value);
    //void SetFrequency(double f);
    //double frequency;
    //Ptr<RandomVariableStream> randomSend;
    std::map<double, std::string> intervalmap;
    std::vector<double> intervalp;
    virtual void ScheduleNextPacket() override;

  protected:
    double rate;

    Time rtt;
    Time MI;

    //utilfunc func;
    double epsilon;
    double theta0;
    double omega0;
    double delta;

    int omegaboundk;


    uint32_t sendnum[2];
    std::unordered_set<uint32_t> sqrMI[2];
    std::queue<Time> delayMI[2];
    std::queue<Time> delayTimestamp[2];

    int8_t probestate;
    int8_t monitorstate;
    EventId probeevent;
    EventId monitorevent;

    bool lastdiffpositive;
    double contk;
};

}  // namespace ndn
}  // namespace ns3
