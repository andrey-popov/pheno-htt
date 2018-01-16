#pragma once

#include <TFile.h>

#include <memory>
#include <queue>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>


class Plugin;


/**
 * \class Processor
 * 
 * Manages a set of plugins and executes them for all events in input files
 * 
 * An object of this class opens input files one by one and feeds them to a set of registered
 * plugins. Plugins are organized into an ordered path and executed consequently for every event.
 * When processing an event, a plugin can reject it, in case of which execution of all subsequent
 * plugins in the path is skipped for this event. A plugin can also notify Processor that there are
 * no more events left in the current input file.
 * 
 * When requested by user, Processor also creates an output ROOT file for each input file.
 * Arbitrary ROOT objects can be created and stored in it.
 */
class Processor
{
public:
    /// Constructor from a collection of paths to input files
    template<typename InputIt>
    Processor(InputIt const &inputFilesBegin, InputIt const &inputFilesEnd);
    
public:
    /**
     * Creates a ROOT object (such as TTree or TH1D) in the output file
     * 
     * A ROOT object is created in the given in-file directory (use "" for the root). Arguments
     * are directly forwarded to the constructor. The returned object is owned by the caller.
     * Objects must be recreated whenever a new input file is opened.
     */
    template<typename T, typename ... Args>
    T *Book(std::string const &inFileDirectory, Args const &... args);
    
    /**
     * Requests automatic creation of output ROOT files
     * 
     * The files will be created in the given directory and named after input files.
     */
    void SetOutput(std::string const outputDir = "");
    
    /**
     * Registers a new plugin
     * 
     * The plugin is added at the end of the path. It is not owned by the instance of Processor.
     * A pointer to Processor is provided to the plugin.
     */
    void RegisterPlugin(Plugin *plugin);
    
    /// Processes input files
    void Run();
    
private:
    /// Opens next input file and, if requested, creates the corresponding output file
    bool OpenNextFile();
    
private:
    /// Queue of input files
    std::queue<std::string> inputFiles;
    
    /// Registered plugins organized in an ordered sequence
    std::vector<Plugin *> path;
    
    /// Current input and output files
    std::unique_ptr<TFile> curInputFile, curOutputFile;
    
    /// Flag showing whether an output ROOT file needs to be created
    bool createOutputFile;
    
    /// Directory in which output ROOT files will be placed
    std::string outputDir;
};


template<typename InputIt>
Processor::Processor(InputIt const &inputFilesBegin, InputIt const &inputFilesEnd):
    createOutputFile(false)
{
    // Save paths to input files
    for (InputIt iter = inputFilesBegin; iter != inputFilesEnd; ++iter)
        inputFiles.push(*iter);
}


template<typename T, typename ... Args>
T *Processor::Book(std::string const &inFileDirectory, Args const &... args)
{
    // Make sure an output file exists
    if (not createOutputFile or not curOutputFile)
    {
        std::ostringstream message;
        message << "Processor::Book: Creation of an output file has not been requested.";
        throw std::runtime_error(message.str());
    }
    
    
    // Change to the given directory. Create it if needed.
    TDirectory *d = curOutputFile->GetDirectory(inFileDirectory.c_str());
    
    if (not d)
        d = curOutputFile->mkdir(inFileDirectory.c_str());
    
    d->cd();
    
    
    // Create the ROOT object in the new current directory
    T *object = new T(args...);
    
    return object;
}
