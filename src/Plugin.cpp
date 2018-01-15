#include <Plugin.hpp>


Plugin::Plugin():
    master(nullptr)
{}


void Plugin::BeginFile(TFile *)
{}


void Plugin::EndFile()
{}


void Plugin::Initialize()
{}


void Plugin::SetMaster(Processor *master_)
{
    master = master_;
}
