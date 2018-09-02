#ifndef __SJC_DEBUG_H__
#define __SJC_DEBUG_H__
#include <iostream>
#include <cassert>
namespace sjc_debug
{
	inline void debugP(size_t& depth, std::string current, std::string destination)
	{
		for(size_t i=0;i!=depth;++i)
			std::cout << "  ";	
		std::cout << "Step into " << destination << " from " << current << std::endl;
		depth++;
	}

	inline void debug0(size_t& depth, std::string content)
	{
		for(size_t i=0;i!=depth;++i)
			std::cout << "  ";	
		std::cout << content << std::endl;
	}

	inline void debugN(size_t& depth, std::string current, std::string destination)
	{
		assert(depth!=0);
		depth--;
		for(size_t i=0;i!=depth;++i)
			std::cout << "  ";	
		std::cout << "Step out of " << destination << " to " << current << std::endl;
	}
}
#endif
