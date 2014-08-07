use output_ops, only: print_err_msg
#define CHCKERR(s) if(ierr.ne.0) then; call print_err_msg(s,rout_name); return; end if
