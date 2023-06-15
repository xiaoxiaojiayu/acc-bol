#ifndef RSTRATEGY
#define RSTRATEGY

#include <limits>
#include <random>

#include "ns3/ndnSIM/NFD/daemon/fw/strategy.hpp"
#include "ns3/ndnSIM/ndn-cxx/util/random.hpp"

namespace nfd {
namespace fw {

using ndn::random::generateWord64;
constexpr uint64_t MAX64 = std::numeric_limits<uint64_t>::max();

class RandomCacheStrategy : public Strategy {
  public:
    explicit RandomCacheStrategy(Forwarder& forwarder,
                          const Name& name = getStrategyName());

    static const Name& getStrategyName();

    void afterReceiveInterest( const FaceEndpoint& ingress, const Interest& interest,
                                            const shared_ptr<pit::Entry>& pitEntry) override;

  private:
    Forwarder& fwd;
    uint64_t boundary2app;
};

}  // namespace fw
}  // namespace nfd
#endif