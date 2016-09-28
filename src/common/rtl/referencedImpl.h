#pragma once

#include "../iface/referencedIface.h"

//#if __cplusplus <= 199711L
    //such a check should work, yet it does not with MS
#if _MSC_VER > 1600L
    #define cpp0x
#else
    #include <Windows.h>
#endif

#ifdef cpp0x
    #include <atomic>
    typedef std::atomic<ULONG> refcounter_t;
#else
    typedef ULONG refcounter_t;
#endif

class CReferenced : public virtual IReferenced
{
    public:
        CReferenced() : mCounter(0) {};
        virtual ~CReferenced() {};

        virtual HRESULT IfaceCalling QueryInterface(/*REFIID */ void*  riid, void ** ppvObj);
        virtual ULONG IfaceCalling AddRef();
        virtual ULONG IfaceCalling Release();

    protected:
        refcounter_t mCounter;
};

class CNotReferenced : public CReferenced
{
    public:
        virtual ULONG IfaceCalling AddRef();
        virtual ULONG IfaceCalling Release();
};
