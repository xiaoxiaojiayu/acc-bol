#include <queue>
#include <unordered_set>

#include "ns3/ndnSIM/apps/ndn-consumer-cbr.hpp"

namespace ns3 {
namespace ndn {

enum PCCCollect_Calc_MODE {
    ENTER_MONITOR_NONE,
    ENTER_MONITOR_RATE_PLUS,
    ENTER_MONITOR_RATE_MINUS,
    ENTER_CLAC_RATE_AND_PROBE
};

enum PCCProbe_MODE {
    ENTER_PROBE_RATE_PAYLOAD_0,
    ENTER_PROBE_RATE_PLUS,
    ENTER_PROBE_RATE_PAYLOAD_1,
    ENTER_PROBE_RATE_MINUS,
    ENTER_PROBE_RATE_PAYLOAD_2
};

typedef struct {
    uint32_t seq;
    int64_t src;
    uint64_t congmark;
    int64_t timestamp;  // relative time point, nanosecond
    int64_t rtt;        // duration, nanosecond
    int64_t rtt_inc_rto;
    bool is_in_sqrMI;
} datainfo_t;

typedef struct {
    double util;
    double lossrate;
} utilret;

typedef double (*utilfunc)(double rate, double drtt_dt, double lrate);

class ConsumerPCC5Stage : public ConsumerCbr {
  public:
    ConsumerPCC5Stage();
    virtual ~ConsumerPCC5Stage();

    static TypeId GetTypeId();

  protected:
    void resettable();

    void updateMaxMinRttFromDatainfo(int index);

    void updateMI();

    // 0: r 1: r-er, 2: r+er
    void PCCProbe(PCCProbe_MODE mode);

    double funcA(double rate, double drtt_dt1, double weight1, double drtt_dt2,
                 double weight2, double markrate);

    double funcB(double rate, double drtt_dt1, double weight1, double drtt_dt2,
                 double weight2, double markrate1,double markrate2);
    double func1(double rate, double mkate, double lrate);

    double funcC(double rate, double drtt_dt1, double weight1, double drtt_dt2,
                 double weight2, double lrate);

    double m(const double &diff);

    utilret calcutil(int i);

    // 0:prepare 1: r-er 2: r+er
    void PCCCollect_Calc(PCCCollect_Calc_MODE mode);

    void ScheduleNextPacket() override;

    void WillSendOutInterest(uint32_t sequenceNumber) override;

    void OnData(shared_ptr<const Data> data) override;

    // negative number: no explicit feedback of src, unknown
    int64_t getSrcofData(shared_ptr<const Data> data);

    void getdelay_inc_rto(Ptr<App> app, uint32_t seqno, Time delay,
                          uint32_t retxCount, int32_t hopCount);

    void getdelay_no_rto(Ptr<App> app, uint32_t seqno, Time delay,
                         int32_t hopCount);

    void SetRandomPrefix(std::string random_prefix);

    void ssrate();

    void StartApplication() override;

  protected:
    double rate;

    Time maxrtt;
    Time minrtt;

    Time rttdiff;
    Time MI;

    // utilfunc func;
    double epsilon;
    double theta0;
    double omega0;
    double delta;
    double v;
    double eta;

    int omegaboundk;
    int marknum = 0;
    int marksum = 0;
    int lossnum = 0;
    double lossrate;
    double randomloss;
    int* queuestart_rec;
    int* queueold_rec;
    int* queuestart_inc;
    int* queueold_inc;
    int srcnum;

    uint32_t sendnum[3];
    std::unordered_set<uint32_t> sqrMI[3];
    std::unordered_map<uint32_t, datainfo_t> datainfo[3];

    int8_t probe_table_index;
    int8_t monitor_table_index;
    EventId probeevent;
    EventId monitorevent;

    bool lastdiffpositive;
    double contk;

    bool randprefix;
    std::map<double, std::string> intervalmap;
    std::vector<double> intervalp;
  
    bool ssfirst;
    int SSflag;
    int decflag;
    int firsttimeset;
    int64_t SStime;
    int64_t dectime;
    double insmarkrate;
    double mrate0;
    std::vector<int> varytest[2];
    int heavycong;
    int receivenum;

    double congestionplus;
    double congestionminus;
     int targetqueue;
    double alpha;
    double beta;
    double gamma;
};

}  // namespace ndn
}  // namespace ns3
