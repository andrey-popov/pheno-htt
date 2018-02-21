/**
 * This program produces ROOT tuples with parton-level and randomly smeared mass of the tt system.
 */

#include <DelphesReaderGen.hpp>
#include <LJetsSelection.hpp>
#include <Processor.hpp>
#include <SmearMttWriter.hpp>

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
    
    DelphesReaderGen reader;
    processor.RegisterPlugin(&reader);
    
    LJetsSelection selection(&reader);
    processor.RegisterPlugin(&selection);
    
    SmearMttWriter writer(&reader, 0.15);
    processor.RegisterPlugin(&writer);
    
    processor.Run();
    
    
    return EXIT_SUCCESS;
}
