#include <queue>
#include <unordered_set>

#include "ns3/ndnSIM/apps/ndn-consumer-cbr.hpp"

namespace ns3 {
namespace ndn {

enum PCCCollect_Calc_MODE {
    ENTER_MONITOR_RATE_PLUS,
    ENTER_MONITOR_RATE_MINUS,
    ENTER_CLAC_RATE_AND_PROBE
};

enum PCCProbe_MODE {
    ENTER_PROBE_RATE,
    ENTER_PROBE_RATE_PLUS,
    ENTER_PROBE_RATE_MINUS
};

typedef struct {
    uint32_t seq;
    int64_t src;
    uint64_t congmark;
    int64_t timestamp;  // relative time point, nanosecond
    int64_t rtt;       // duration, nanosecond
    bool is_in_sqrMI;
} datainfo_t;

typedef struct {
    double util;
    double avgrtt;
} utilret;

typedef double (*utilfunc)(double rate, double drtt_dt, double lrate);

class ConsumerPCCMY : public ConsumerCbr {
  public:
    ConsumerPCCMY();
    virtual ~ConsumerPCCMY();

    static TypeId GetTypeId();

    void setrate(double ir);

  protected:
    void changeprefix();

    void resettable();

    void updateMI();

    // 0: r 1: r-er, 2: r+er
    void PCCProbe(uint8_t mode);

    double funcA(double rate, double drtt_dt1, double weight1, double drtt_dt2,
                double weight2, double markrate);
    
    double funcB(double rate, double drtt_dt1, double weight1, double drtt_dt2,
                double weight2, double markrate);
    double func1(double rate, double mkate, double lrate);
    
    double funcC(double rate, double drtt_dt1, double weight1, double drtt_dt2,
                double weight2, double lrate);

    double m(const double &diff);

    utilret calcutil(int i);

    // 0:prepare 1: r-er 2: r+er
    void PCCCollect_Calc(uint8_t mode);

    void ScheduleNextPacket() override;

    void WillSendOutInterest(uint32_t sequenceNumber) override;

    void OnData(shared_ptr<const Data> data) override;

    // negative number: no explicit feedback of src, unknown
    int64_t getSrcofData(shared_ptr<const Data> data);

    void getdelay(Ptr<App> app, uint32_t seqno, Time delay, uint32_t retxCount,
                  int32_t hopCount);

    void SetRandomPrefix(std::string random_prefix);
   
    void ssrate ();

    void StartApplication() override;

    virtual void OnTimeout(uint32_t sequenceNum) override;
    //virtual void OnTimeout(uint32_t sequenceNum) override;

  protected:
    double rate;

    Time minrtt;
    Time maxrtt;
    Time rttdiff;
    Time rtt;
    Time MI;
    int SSflag;
    int decflag;
    int fincreaseflag;
    int64_t SStime;
    int64_t dectime;
    double insmarkrate;
    double mrate0;
    std::vector<int> varytest[2];
    int heavycong;
    int receivenum;
    double timeoutnum_minus;
    double timeoutnum_plus;
    double timeoutnum;
    int stateflag;

    // utilfunc func;
    double epsilon;
    double theta0;
    double omega0;
    double delta;

    int omegaboundk;
    int marknum = 0;
    int marksum = 0;
    int lossnum = 0;
   
    uint32_t sendnum[3];
    std::unordered_set<uint32_t> sqrMI[3];
    std::unordered_map<uint32_t, datainfo_t> datainfo[3];

    int8_t probestate;
    int8_t monitorstate;
    EventId probeevent;
    EventId monitorevent;

    bool lastdiffpositive;
    double contk;
    int srcnum;

    bool randprefix;
    std::map<double, std::string> intervalmap;
    std::vector<double> intervalp;
};

}  // namespace ndn
}  // namespace ns3
