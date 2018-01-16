#include <DelphesReader.hpp>
#include <LJetsSelection.hpp>
#include <Processor.hpp>

#include <iostream>


using std::cout;


int main(int argc, char **argv)
{
    if (argc != 2)
    {
        std::cerr << "Usage: htt-tuples INPUT_FILE_MASK\n";
        return EXIT_FAILURE;
    }
    
    
    Processor processor(argv[1]);
    
    DelphesReader reader;
    processor.RegisterPlugin(&reader);
    
    LJetsSelection selection(&reader);
    processor.RegisterPlugin(&selection);
    
    processor.Run();
    
    
    return EXIT_SUCCESS;
}
