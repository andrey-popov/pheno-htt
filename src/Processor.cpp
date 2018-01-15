#include <Processor.hpp>

#include <Plugin.hpp>

#include <boost/filesystem.hpp>

#include <algorithm>
#include <iostream>


Plugin *Processor::GetPlugin(std::string const &name, bool checkState)
{
    auto const res = std::find_if(path.begin(), path.end(),
      [&name](PluginInPath const &p){return (p.name == name);});
    
    if (res == path.end())
    {
        std::ostringstream message;
        message << "Processor::GetPlugin: No plugin with name \"" << name <<
          "\" has been registered";
        throw std::runtime_error(message.str());
    }
    
    if (checkState and not res->hasRun)
    {
        std::ostringstream message;
        message << "Processor::GetPlugin: Requested plugin with name \"" << name <<
          "\" has not been run for the current event yet.";
        throw std::runtime_error(message.str());
    }
    
    return res->plugin.get();
}


void Processor::SetOutput(std::string const outputDir_)
{
    createOutputFile = true;
    outputDir = outputDir_;
}


void Processor::RegisterPlugin(std::unique_ptr<Plugin> &&plugin, std::string const &name)
{
    // Make sure no plugin with the same name has been registered
    auto res = std::find_if(path.begin(), path.end(),
      [&name](PluginInPath const &p){return (p.name == name);});
    
    if (res != path.end())
    {
        std::ostringstream message;
        message << "Processor::RegisterPlugin: Cannot register a plugin with name \"" << name <<
          "\" because another plugin with the same name has already been registered.";
        throw std::runtime_error(message.str());
    }
    
    
    // Register the plugin, also giving it a link to Processor
    plugin->SetMaster(this);
    path.emplace_back(PluginInPath{std::move(plugin), name, false});
}


void Processor::Run()
{
    // Create directory for output files if needed
    if (createOutputFile)
        boost::filesystem::create_directories(outputDir);
    
    
    // Initialization for plugins
    for (auto &p: path)
        p.plugin->Initialize();
    
    
    // Loop over input files
    while (OpenNextFile())
    {
        // Notify all plugins that a new input file has been started
        for (auto &p: path)
            p.plugin->BeginFile(curInputFile.get());
        
        
        // Loop over events in the current input file
        bool noEvents = false;
        
        while (not noEvents)
        {
            // Reset status flags for all plugins
            for (auto &p: path)
                p.hasRun = false;
            
            // Process current event through all plugins
            for (auto &p: path)
            {
                Plugin::EventOutcome const res = p.plugin->ProcessEventToOutcome();
                p.hasRun = true;
                
                if (res == Plugin::EventOutcome::NoEvents)
                {
                    noEvents = true;
                    break;
                }
            }
        }
        
        
        // Notify all plugins that the current input file is about to be closed
        for (auto pIt = path.rbegin(); pIt != path.rend(); ++pIt)
            pIt->plugin->EndFile();
    }
}


bool Processor::OpenNextFile()
{
    // Close previous input and output files
    if (curInputFile)
        curInputFile->Close();
    
    if (curOutputFile)
    {
        curOutputFile->Write();
        curOutputFile->Close();
    }
    
    
    // Check if there are input files left
    if (inputFiles.empty())
        return false;
    
    
    // Open the next input file
    std::string const inputFileName(inputFiles.front());
    curInputFile.reset(TFile::Open(inputFileName.c_str()));
    
    if (not curInputFile or curInputFile->IsZombie())
    {
        std::ostringstream message;
        message << "Processor::OpenNextFile: Failed to open file \"" << inputFileName <<
          "\" for reading.";
        throw std::runtime_error(message.str());
    }
    
    std::cout << "Processing file \"" << inputFileName << "\"..." << std::endl;
    inputFiles.pop();
    
    
    // Create output file
    if (createOutputFile)
    {
        auto const inputBaseName = boost::filesystem::path(inputFileName).filename();
        auto const outFileName = (boost::filesystem::path(outputDir) / inputBaseName).string();
        
        curOutputFile.reset(TFile::Open(outFileName.c_str(), "create"));
        
        if (not curOutputFile or curOutputFile->IsZombie())
        {
            std::ostringstream message;
            message << "Failed to open file \"" << outFileName << "\" for writing.";
            
            if (boost::filesystem::exists(outFileName))
                message << " The file already exists.";
            
            throw std::runtime_error(message.str());
        }
    }
    
    
    return true;
}
