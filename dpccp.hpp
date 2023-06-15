#ifndef DPCCP
#define DPCCP
#include <queue>
#include <unordered_set>

#include "ns3/ndnSIM/apps/ndn-consumer-cbr.hpp"

namespace ns3 {
namespace ndn {

constexpr double MIN_RATE = 2;

class ConsumerDPCCP : public ConsumerCbr {
  private:
    enum class States;
    enum class BFStates;
    enum class FCStates;

    typedef struct {
        uint32_t seq;
        int64_t respondid;
        int64_t rtt;
        int64_t timestamp;
    } datainfo;

    typedef struct {
        double backlogged;
        double receive_rate;
    } bcrettype;

  public:
    ConsumerDPCCP();
    virtual ~ConsumerDPCCP();

    static TypeId GetTypeId();

  protected:
    void ProbeRtt();

    void BFEntry();
    void BFDoing(BFStates);

    void FCEntry();
    void FCDoing(FCStates);

    void QDEntry();
    void QDDoing();

    void OnData(shared_ptr<const Data> data) override;
    void ScheduleNextPacket() override;

    void getdelay(Ptr<App> app, uint32_t seqno, Time delay, 
                  int32_t hopCount);
     void SetRandomPrefix(std::string random_prefix);

    // void WillSendOutInterest(uint32_t sequenceNumber) override;

    // void ScheduleNextPacket() override;

    int64_t getSrcofData(shared_ptr<const Data> data);

    bcrettype calc_backlogged_by_datainfomap();

       void StartApplication() override;

  private:
    friend std::ostream& operator<<(std::ostream& os, States s);
    friend std::ostream& operator<<(std::ostream& os, BFStates s);
    friend std::ostream& operator<<(std::ostream& os, FCStates s);

  protected:
    States curstate;

    // unit: pkts/s
    double rate;

    // rate realated parameters
    double beta;
    double zeta;
    int SSflag;
    // threshold
    double psi;
    double alpha_r;
    int64_t rho;
    double eta;
    bool randprefix;
    int64_t settime0;
    int64_t settime1;
    std::map<double, std::string> intervalmap;
    std::vector<double> intervalp;

    // states-transfer related time
    Time baseRtt[5];
    Time baseRtt0;
    Time baseRtt0_step;
    Time baseRtt1;
    Time baseRtt1_step;
    Time MI;
    Time maxrtt;
    int srcnum;
    int rttupdateflag;
    std::unordered_map<uint32_t, datainfo> datainfomap;
    // std::unordered_set<int64_t> rtt_by_src[2];

    double fc_delta_r[3];

    // bool isBF2FC;
    // double receive_rate;

    EventId nextevent;
};

}  // namespace ndn
}  // namespace ns3

#endif
