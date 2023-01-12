use parameters;
use QGN_Module;

proc main() {

  Initialize();

  /*
  for i in 1..Nt {

    call TimeStep(q_hat,t,dt)
    if( mod(i,out_freq)==0 ) print *, 'Time since inception = ',real(t/86400._dp),'days; time step size = ',real(dt),'s'
    if( mod(i,diag_freq)==0 ) then
        print *, ' Diagnostics, RecNo = ', i/diag_freq
        call WritePUV(q_hat,i/diag_freq) ! Writes psi to pSeries.dat, also u and v
        call writeDiagnostics(i/diag_freq) ! Writes snapshots of KE budget terms
    end if
    if( mod(i,chkpt_freq)==0 ) then
        call WriteQP(q_hat,1) ! Writes q and psi to q.dat and p.dat.
                              ! Good for restarts.
    end if
end do
call Cleanup

end program main
     */
}
