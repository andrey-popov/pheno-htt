/**
 * This program evaluates performance of tt reconstruction.
 */

#include <DelphesReader.hpp>
#include <LJetsSelection.hpp>
#include <Processor.hpp>
#include <TTReco.hpp>
#include <TTRecoPerf.hpp>

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
    
    TTReco ttReco(&reader, &selection, "data/tt-reco.root");
    processor.RegisterPlugin(&ttReco);
    
    TTRecoPerf perf(&reader, &ttReco, "tt-reco-performance.root");
    processor.RegisterPlugin(&perf);
    
    processor.Run();
    
    cout << '\n';
    perf.PrintCounts();
    
    
    return EXIT_SUCCESS;
}
