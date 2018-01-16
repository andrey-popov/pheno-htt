/**
 * This program produces histograms that are needed for evaluation of the likelihood function used
 * in reconstruction of tt -> l+jets decays.
 */

#include <DelphesReader.hpp>
#include <LJetsSelection.hpp>
#include <Processor.hpp>
#include <TTRecoInputs.hpp>

#include <iostream>


using std::cout;


int main(int argc, char **argv)
{
    if (argc < 2)
    {
        std::cerr << "Usage: htt-tuples INPUT_FILE_MASK\n";
        return EXIT_FAILURE;
    }
    
    
    Processor processor(argv + 1, argv + argc);
    
    DelphesReader reader(DelphesReader::LHE_PARTICLES);
    processor.RegisterPlugin(&reader);
    
    LJetsSelection selection(&reader);
    processor.RegisterPlugin(&selection);
    
    TTRecoInputs recoBuilder(&reader, &selection, "tt-reco.root");
    processor.RegisterPlugin(&recoBuilder);
    
    processor.Run();
    
    cout << '\n';
    recoBuilder.PrintCounts();
    
    
    return EXIT_SUCCESS;
}
