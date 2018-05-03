/**
 * This program produces ROOT histograms with randomly smeared mass of the tt system.
 */

#include <DelphesReaderGen.hpp>
#include <LJetsLHEFilter.hpp>
#include <LJetsSelection.hpp>
#include <Processor.hpp>
#include <SystMttHists.hpp>

#include <boost/program_options.hpp>

#include <iostream>
#include <vector>


int main(int argc, char **argv)
{
    namespace po = boost::program_options;
    
    // Parse arguments
    po::options_description options("Allowed options");
    options.add_options()
      ("help,h", "Prints help message")
      ("inputFiles", po::value<std::vector<std::string>>(), "Input files")
      ("resolution,r", po::value<double>()->default_value(0.2), "Relative mtt resolution");
    
    po::positional_options_description positionalOptions;
    positionalOptions.add("inputFiles", -1);
    
    po::variables_map optionsMap;
    
    po::store(
      po::command_line_parser(argc, argv).options(options).positional(positionalOptions).run(),
      optionsMap);
    po::notify(optionsMap);
    
    if (optionsMap.count("help"))
    {
        std::cerr << "Produces histograms with smeared mass of tt system.\n";
        std::cerr << "Usage: mtt-hists inputFiles [options]\n";
        std::cerr << options << std::endl;
        return EXIT_FAILURE;
    }
    
    if (not optionsMap.count("inputFiles"))
    {
        std::cerr << "Usage: mtt-hists inputFiles [options]\n";
        std::cerr << options << std::endl;
        return EXIT_FAILURE;
    }
    
    std::vector<std::string> inputFiles(optionsMap["inputFiles"].as<std::vector<std::string>>());
    
    
    Processor processor(inputFiles.begin(), inputFiles.end());
    processor.SetOutput("output");
    
    DelphesReaderGen reader;
    reader.SetReadLHEWeights();
    processor.RegisterPlugin(&reader);
    
    LJetsLHEFilter lheFilter(&reader);
    processor.RegisterPlugin(&lheFilter);
    
    LJetsSelection selection(&reader);
    processor.RegisterPlugin(&selection);
    
    
    std::vector<double> binning{350, 368, 388, 408, 430, 452, 476, 501, 528, 556, 585, 616, 648, 682, 718, 756, 796, 838, 882, 928, 977, 1029, 1083, 1140, 1200};
    
    SystMttHists writer(&reader, binning, optionsMap["resolution"].as<double>());
    processor.RegisterPlugin(&writer);
    
    processor.Run();
    
    
    return EXIT_SUCCESS;
}
