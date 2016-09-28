#pragma once

#include <vector>

#include "../../common/iface/ApproxIface.h"
#include "..\..\common\rtl\hresult.h"
#include "..\..\common\rtl\referencedImpl.h"


#pragma warning( push )
#pragma warning( disable : 4250 ) // C4250 - 'class1' : inherits 'class2::member' via dominance

class CGlucoseLevels : public IGlucoseLevels, public virtual CReferenced{
protected:
	std::vector<TGlucoseLevel> mLevels;
public:	
	virtual ~CGlucoseLevels() {};
		//dctor has to be virtual, even if it is empty, due to the inheritance by dominance

	HRESULT IfaceCalling GetLevels(TGlucoseLevel** levels);
	HRESULT IfaceCalling GetLevelsCount(size_t* count);
	HRESULT IfaceCalling SetLevelsCount(size_t count);	
	HRESULT IfaceCalling GetBounds(TGlucoseLevelBounds *bounds);	
};

#pragma warning( pop )