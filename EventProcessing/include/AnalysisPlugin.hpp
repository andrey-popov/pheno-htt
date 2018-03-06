#pragma once

#include <Plugin.hpp>



/**
 * \class AnalysisPlugin
 * 
 * A special version of Plugin that cannot terminate processing of the current input file
 * 
 * To be used for plugins that do not access input files directly. Processing of events must be
 * done by implementing pure virtual method ProcessEvent.
 */
class AnalysisPlugin: public Plugin
{
public:
    virtual Plugin::EventOutcome ProcessEventToOutcome() override final;
    
private:
    /**
     * Processes current event
     * 
     * Returned boolean is interpreted as a request to keep or reject current event.
     */
    virtual bool ProcessEvent() = 0;
};
