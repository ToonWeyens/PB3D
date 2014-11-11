use message_ops, only: print_err_msg
#define CHCKERR(s) if(ierr.ne.0) then; if(ierr.ne.66) call print_err_msg(s,rout_name); return; end if
#define CHCKSTT if(istat.ne.0) then; return; end if
