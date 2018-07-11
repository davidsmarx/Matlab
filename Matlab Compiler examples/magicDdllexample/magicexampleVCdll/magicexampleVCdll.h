// The following ifdef block is the standard way of creating macros which make exporting 
// from a DLL simpler. All files within this DLL are compiled with the MAGICEXAMPLEVCDLL_EXPORTS
// symbol defined on the command line. this symbol should not be defined on any project
// that uses this DLL. This way any other project whose source files include this file see 
// MAGICEXAMPLEVCDLL_API functions as being imported from a DLL, whereas this DLL sees symbols
// defined with this macro as being exported.
#ifdef MAGICEXAMPLEVCDLL_EXPORTS
#define MAGICEXAMPLEVCDLL_API __declspec(dllexport)
#else
#define MAGICEXAMPLEVCDLL_API __declspec(dllimport)
#endif

// This class is exported from the magicexampleVCdll.dll
class MAGICEXAMPLEVCDLL_API CmagicexampleVCdll {
public:
	CmagicexampleVCdll(void);
	// TODO: add your methods here.
};

extern MAGICEXAMPLEVCDLL_API int nmagicexampleVCdll;

MAGICEXAMPLEVCDLL_API int fnmagicexampleVCdll(void);
