#include "referencedImpl.h"

#include "cfixes.h"

#include <limits>

HRESULT CReferenced::QueryInterface(/*REFIID */ void*  riid, void ** ppvObj)
{
    unused(ppvObj);
    unused(riid);
    return E_NOTIMPL;
}

ULONG IfaceCalling CReferenced::AddRef()
{
#ifdef cpp0x
    return mCounter.fetch_add(1) + 1;
#else
    return InterlockedIncrement(&mCounter);
#endif
}

ULONG IfaceCalling CReferenced::Release()
{
#ifdef cpp0x
    ULONG rc = mCounter.fetch_sub(1) - 1;   //fetch_sub returns the old value!
#else
    ULONG rc = InterlockedDecrement(&mCounter);
#endif
    if (rc == 0) delete this;

    return rc;
}

#undef max

ULONG IfaceCalling CNotReferenced::AddRef()
{
    return std::numeric_limits<ULONG>::max();
}

ULONG IfaceCalling CNotReferenced::Release()
{
    return std::numeric_limits<ULONG>::max();
}
