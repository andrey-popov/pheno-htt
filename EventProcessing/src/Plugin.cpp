#include <Plugin.hpp>


Plugin::Plugin():
    processor(nullptr)
{}


void Plugin::BeginFile(TFile *)
{}


void Plugin::EndFile()
{}


void Plugin::SetProcessor(Processor *processor_)
{
    processor = processor_;
}
