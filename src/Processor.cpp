#include <Processor.hpp>

#include <Plugin.hpp>

#include <boost/filesystem.hpp>

#include <algorithm>
#include <iostream>


void Processor::SetOutput(std::string const outputDir_)
{
    createOutputFile = true;
    outputDir = outputDir_;
}


void Processor::RegisterPlugin(Plugin *plugin)
{
    plugin->SetProcessor(this);
    path.emplace_back(plugin);
}


void Processor::Run()
{
    // Create directory for output files if needed
    if (createOutputFile)
        boost::filesystem::create_directories(outputDir);
    
    
    // Loop over input files
    while (OpenNextFile())
    {
        // Notify all plugins that a new input file has been started
        for (auto &p: path)
            p->BeginFile(curInputFile.get());
        
        
        // Loop over events in the current input file
        bool noEvents = false;
        
        while (not noEvents)
        {
            // Process current event through all plugins
            for (auto &p: path)
            {
                Plugin::EventOutcome const res = p->ProcessEventToOutcome();
                
                if (res == Plugin::EventOutcome::NoEvents)
                {
                    noEvents = true;
                    break;
                }
            }
        }
        
        
        // Notify all plugins that the current input file is about to be closed
        for (auto pIt = path.rbegin(); pIt != path.rend(); ++pIt)
            (*pIt)->EndFile();
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
