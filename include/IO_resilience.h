#if lrIO
use num_vars, only: mnt => max_nr_tries_IO, net => nr_extra_tries_IO, rIOd
#define rIO(a,b) do rIOd = 1,mnt; a; if (b.eq.0) exit; call sleep(1); net = net+1; end do
#define rIO2(a,b) do rIOd = 1,mnt; a; if (b.eq.0) exit; call sleep(1); end do
#else
#define rIO(a,b) a
#define rIO2(a,b) a
#endif
