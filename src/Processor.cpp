#include <Processor.hpp>

#include <Plugin.hpp>

#include <boost/filesystem.hpp>

#include <algorithm>
#include <iostream>
#include <regex>


Processor::Processor(std::string const &fileMask_):
    createOutputFile(false)
{
    namespace fs = boost::filesystem;
    
    // Split the mask into directory and filename parts and make sure that the former one do not
    //contain wildcards
    fs::path const fileMask(fileMask_);
    fs::path const directory(fileMask.parent_path());
    
    if (directory.string().find_first_of("*?") != std::string::npos)
    {
        std::ostringstream message;
        message << "Processor::Processor: Directory part of pattern \"" << directory.string() <<
          "\" contains wildcards, which is not supported.";
        throw std::runtime_error(message.str());
    }
    
    
    // Make sure the directory exists
    if (not fs::exists(directory) or not fs::is_directory(directory))
    {
        std::ostringstream message;
        message << "Processor::Processor: Directory \"" << directory.string() <<
          "\" does not exist.";
        throw std::runtime_error(message.str());
    }
    
    
    // Convert the filename pattern into a regular expression. In order to do it, escape all
    //special characters except for '*' and '?' and prepend these two with dots.
    std::string filenameMask(fileMask.filename().string());
    
    std::regex escapeRegex(R"(\.|\^|\$|\||\(|\)|\[|\]|\{|\}|\+|\\)");
    filenameMask = std::regex_replace(filenameMask, escapeRegex, "\\$&");
    
    std::regex prependRegex(R"(\?|\*)");
    filenameMask = std::regex_replace(filenameMask, prependRegex, ".$&");
    
    std::regex filenameRegex(filenameMask);
    
    
    // Check all files in the directory against the constructed regular expression
    for (fs::directory_iterator dirIt(directory); dirIt != fs::directory_iterator(); ++dirIt)
    {
        if (not fs::is_regular_file(dirIt->status()))
            continue;
        
        if (std::regex_match(dirIt->path().filename().string(), filenameRegex))
            inputFiles.push(dirIt->path().string());
    }
    
    if (inputFiles.empty())
    {
        std::ostringstream message;
        message << "Processor::Processor: Found no file matching mask \"" << fileMask_ << "\".";
        throw std::runtime_error(message.str());
    }
}


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
