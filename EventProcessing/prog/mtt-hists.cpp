/**
 * This program produces ROOT histograms with randomly smeared mass of the tt system.
 */

#include <DelphesReaderGen.hpp>
#include <LJetsLHEFilter.hpp>
#include <LJetsSelection.hpp>
#include <Processor.hpp>
#include <SystMttHists.hpp>

#include <iostream>
#include <vector>


int main(int argc, char **argv)
{
    if (argc < 2)
    {
        std::cerr << "Usage: htt-tuples INPUT_FILE_MASK\n";
        return EXIT_FAILURE;
    }
    
    
    Processor processor(argv + 1, argv + argc);
    processor.SetOutput("output");
    
    DelphesReaderGen reader;
    reader.SetReadLHEWeights();
    processor.RegisterPlugin(&reader);
    
    LJetsLHEFilter lheFilter(&reader);
    processor.RegisterPlugin(&lheFilter);
    
    LJetsSelection selection(&reader);
    processor.RegisterPlugin(&selection);
    
    
    std::vector<double> binning{350, 368, 388, 408, 430, 452, 476, 501, 528, 556, 585, 616, 648, 682, 718, 756, 796, 838, 882, 928, 977, 1029, 1083, 1140, 1200};
    
    SystMttHists writer(&reader, binning, 0.2);
    processor.RegisterPlugin(&writer);
    
    processor.Run();
    
    
    return EXIT_SUCCESS;
}
