#if ( lwith_intel && !lwith_gnu)
#define rp real
#define ip aimag
#else
#define rp realpart
#define ip imagpart
#endif
