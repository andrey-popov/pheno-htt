#pragma once


class Processor;
class TFile;


/**
 * \class Plugin
 * 
 * Abstract base class to describe plugins
 * 
 * A derived class must implement a method to process an event. It can also optionally implement
 * methods to be called whenever a new input file is opened or closed.
 */
class Plugin
{
public:
    /// Possible outcomes of processing of a single event
    enum class EventOutcome
    {
        Ok,        ///< Everything is fine
        Rejected,  ///< Event is rejected
        NoEvents   ///< There are no more events in the input file
    };
    
public:
    Plugin();
    
    virtual ~Plugin() = default;
    
public:
    /**
     * Notifies the plugin that a new input file has been opened
     * 
     * A non-owning pointer to the file is given as the argument. Default implementation is
     * trivial.
     */
    virtual void BeginFile(TFile *inputFile);
    
    /**
     * Notifies the plugin that the current input file is about to be closed
     * 
     * Default implementation is trivial.
     */
    virtual void EndFile();
    
    /// Requests processing of the current event
    virtual EventOutcome ProcessEventToOutcome() = 0;
    
    /// Sets non-owning pointer to an instance of Processor to which this plugin is attached
    void SetProcessor(Processor *processor);
    
private:
    /**
     * Non-owning pointer to Processor to which this plugin is attached
     * 
     * Needed to allow communication with other plugins.
     */
    Processor *processor;
};
