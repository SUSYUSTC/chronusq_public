#include <singleslater.hpp>
#include <iostream>
namespace ChronusQ
{
	void SingleSlaterBase::initpcm(std::shared_ptr<PCMBase> pcmbase)
	{
		this->pcm=pcmbase;
	}
}
