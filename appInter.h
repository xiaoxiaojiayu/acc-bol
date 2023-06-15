#ifndef APP_INTER_H
#define APP_INTER_H
#include "ns3/ndnSIM/model/ndn-common.hpp"

#include "ns3/ndnSIM/apps/ndn-app.hpp"
#include "ns3/ndnSIM/model/ndn-common.hpp"
#include "ns3/ndnSIM/model/ndn-l3-protocol.hpp"
#include "ns3/random-variable-stream.h"



#include "ns3/ndnSIM/NFD/daemon/fw/forwarder.hpp"
#include "ns3/ndnSIM/NFD/daemon/face/internal-face.hpp"
#include "ns3/ndnSIM/NFD/daemon/face/internal-transport.hpp"
#include "ns3/ndnSIM/NFD/daemon/mgmt/fib-manager.hpp"
#include "ns3/ndnSIM/NFD/daemon/mgmt/face-manager.hpp"
#include "ns3/ndnSIM/NFD/daemon/mgmt/strategy-choice-manager.hpp"
#include "ns3/ndnSIM/NFD/daemon/mgmt/cs-manager.hpp"
#include "ns3/ndnSIM/NFD/daemon/mgmt/forwarder-status-manager.hpp"
#include "ns3/ndnSIM/NFD/daemon/mgmt/tables-config-section.hpp"
#include "ns3/ndnSIM/NFD/daemon/mgmt/command-authenticator.hpp"

#include "ns3/ndnSIM/NFD/daemon/rib/service.hpp"

#include "ns3/ndnSIM/NFD/daemon/face/null-face.hpp"
#include "ns3/ndnSIM/NFD/daemon/face/internal-face.hpp"

#include "ns3/ndnSIM/NFD/daemon/common/global.hpp"
#include "ns3/ndnSIM/NFD/daemon/common/config-file.hpp"

#include "ns3/nstime.h"
#include "ns3/ptr.h"
namespace ns3 {
    namespace ndn {

        class AppInter : public App {
        public:
            static TypeId GetTypeId(void);
            AppInter();

            virtual void OnInterest(shared_ptr<const Interest> interest);

        protected:
            virtual void StartApplication();
            virtual void StopApplication();

        private:
            Name prefix;
            Ptr<UniformRandomVariable> rand;
            Ptr<L3Protocol> l3protocol;
            shared_ptr<nfd::Forwarder> forwarder;
        };
    }
}
#endif