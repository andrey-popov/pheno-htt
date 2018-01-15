#include <AnalysisPlugin.hpp>


Plugin::EventOutcome AnalysisPlugin::ProcessEventToOutcome()
{
    bool const result = ProcessEvent();
    
    if (result)
        return Plugin::EventOutcome::Ok;
    else
        return Plugin::EventOutcome::Rejected;
}
