program main
    use variables
    
    implicit none
    real(4)::startt,endt

    startt=secnds(0.0)

    ep = 'EPD' !EPD ou EPT
    numberOfThreads=8
    
    call InputData()
    call Kernel()
    call OutputData()
    !call OutputDataForGfortran()

    endt=secnds(startt)

    !call Cplusplus()

    write(*,'(x,a)',advance='no')'TEMPO DE EXECUCAO: '
    write(*,'(f9.5)',advance='no')endt
    write(*,'(a)')' SEGUNDOS.'
    write(*,*)' '



    pause










end program main