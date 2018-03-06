#include <DelphesReader.hpp>
#include <LJetsSelection.hpp>
#include <Processor.hpp>
#include <TTReco.hpp>
#include <VarWriter.hpp>

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
    processor.SetOutput("output");
    
    DelphesReader reader;
    processor.RegisterPlugin(&reader);
    
    LJetsSelection selection(&reader);
    processor.RegisterPlugin(&selection);
    
    TTReco ttReco(&reader, &selection, "data/tt-reco.root");
    processor.RegisterPlugin(&ttReco);
    
    VarWriter writer(&reader, &ttReco);
    processor.RegisterPlugin(&writer);
    
    processor.Run();
    
    
    return EXIT_SUCCESS;
}
