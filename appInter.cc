#include "ns3/ndnSIM/model/ndn-l3-protocol.hpp"
#include "ns3/ndnSIM/model/ndn-app-link-service.hpp"
#include "ns3/ndnSIM/model/null-transport.hpp"
#include "ns3/ndnSIM/helper/ndn-fib-helper.hpp"
#include "ns3/string.h"
#include "appInter.h"

NS_LOG_COMPONENT_DEFINE("ndn.AppInter");

namespace ns3 {
    namespace ndn {

        NS_OBJECT_ENSURE_REGISTERED(AppInter);

        TypeId AppInter::GetTypeId() {
            static TypeId tid =
                TypeId("ns3::ndn::AppInter")
                .SetGroupName("Ndn")
                .SetParent<App>()
                .AddConstructor<AppInter>()
                .AddAttribute("Prefix", "", StringValue("/"),
                    MakeNameAccessor(&AppInter::prefix), MakeNameChecker());
            return tid;
        }

        AppInter::AppInter()
            : rand(CreateObject<UniformRandomVariable>()) {}

        void AppInter::StartApplication() {
            NS_LOG_FUNCTION_NOARGS();
            App::StartApplication();
            l3protocol = GetNode()->GetObject<L3Protocol>();
            forwarder = l3protocol->getForwarder();
            FibHelper::AddRoute(GetNode(), prefix, m_face, 1);
            NS_LOG_DEBUG("l3:" << l3protocol << "\tfwd:" << forwarder << "\tname:" << prefix);
        }

        void AppInter::StopApplication() {
            NS_LOG_FUNCTION_NOARGS();
            App::StopApplication();
        }

        void AppInter::OnInterest(shared_ptr<const Interest> interest) {
            App::OnInterest(interest);
            NS_LOG_DEBUG("Interest seq:" << interest->getName().get(-1).toSequenceNumber()
                << " PitSize:" << forwarder->getPit().size());

            auto nonconstint = const_cast<Interest*>(&(*interest));
            nonconstint->setNonce(rand->GetValue(0, std::numeric_limits<uint32_t>::max()));

            m_transmittedInterests(interest, this, m_face);
            m_appLink->onReceiveInterest(*interest);
        }
    }
}
