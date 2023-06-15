#include <queue>
#include <unordered_set>

#include "ns3/ndnSIM/apps/ndn-consumer-cbr.hpp"

namespace ns3 {
namespace ndn {


class ConsumerSINR : public ConsumerCbr {
  public:
    ConsumerSINR();
    virtual ~ConsumerSINR();

    static TypeId GetTypeId();

    void setrate(double ir);

  protected:

    void resettable();
    void ScheduleNextPacket() override;
    // 0: r 1: r-er, 2: r+er
    void OnData(shared_ptr<const Data> data) override;
    virtual void OnTimeout(uint32_t sequenceNum) override;

    // negative number: no explicit feedback of src, unknown
   
    void getdelay(Ptr<App> app, uint32_t seqno, Time delay, uint32_t retxCount,
                  int32_t hopCount);
    void SetRandomPrefix(std::string random_prefix);
   
    void ssrate ();
    void fsupdate();

  protected:
    double rate;
    int SSflag;
    double SStime;
    double outtime;
    double outnum;
    double dectime;

    double fsmax;
    double fsmin;
    double fsnew;
    double fsold;
    double fsstep;
    double Terror;
    double Error_IDG;
    //double IDG;
    double IDG_new;
    double IDG_old;
    double IIG;
    int Nidg;
    int cur_num;
    double alpha;
    double beta;
    int beforestart;

    Time MI;
    double rttnew;
    double rttold;
    double rttm;
    double varnew;
    double varold;
    double RTOnew;
    int firststartflag;

    std::vector<int64_t> data_time;
    int64_t data_interval;


    bool randprefix;
    std::map<double, std::string> intervalmap;
    std::vector<double> intervalp;
};

}  // namespace ndn
}  // namespace ns3
