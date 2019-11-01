
subroutine InputData()
    use variables

    implicit none

    integer::iostatus,cont,cont1,cont2,i=0,ih=0
    real(8)::value
    character(6)::teste1, teste2
    real(4)::start, end

    write(*,'(x,a)')'INICIANDO LEITURA DE DADOS.'

    start=secnds(0.0)

    !MPLIST
    open(11,file='MPLIST.lis',status='unknown')
    iostatus=0
    cont=0
    
    do while(iostatus==0)
        read(11,*,iostat=iostatus)teste1
        if(teste1=='EX   ') cont=cont+1     
    end do
    close(11)

    allocate(materials(cont,3))
    materials=0.0
    open(11,file='MPLIST.lis',status='unknown')
    read(11,*,iostat=iostatus)
    read(11,*,iostat=iostatus)
    read(11,*,iostat=iostatus)
    
    iostatus=0 
    cont=0

    do while(iostatus==0)
        read(11,*,iostat=iostatus)teste1
    
        if(teste1=='MATERI') then
            cont=cont+1
            read(11,*,iostat=iostatus)teste1,teste2,materials(cont,1)
            read(11,*,iostat=iostatus)teste1,teste2,materials(cont,2)
            
        end if

    end do 
    close(11)

    !RLIST
    iostatus=0
    open(22,file='RLIST.lis',status='unknown')
    read(22,*,iostat=iostatus)
    read(22,*,iostat=iostatus)
    read(22,*,iostat=iostatus)
    do i=1,size(materials,dim=1)
        read(22,*,iostat=iostatus)
        read(22,*,iostat=iostatus)materials(i,3)
        read(22,*,iostat=iostatus)
    end do
    close(22)

    !NLIST
    cont=0
    iostatus=0
    open(33,file='NLIST.lis',status='unknown')
    do while(iostatus==0)
        read(33,*,iostat=iostatus)teste1      
    end do
    close(33)

    read(teste1(1:6),*)cont
    allocate(coord(cont,2))
    coord=0.0
    open(33,file='NLIST.lis',status='unknown')
    read(33,*,iostat=iostatus)
    read(33,*,iostat=iostatus)
    read(33,*,iostat=iostatus)
    read(33,*,iostat=iostatus)
    cont1=cont/20
    cont2=MOD(cont,20)
    cont=0
    do i=1,cont1
        read(33,*,iostat=iostatus)
        do ih=1,20
            cont=cont+1
            read(33,*,iostat=iostatus)teste1,coord(cont,1),coord(cont,2)
        end do !ih
        read(33,*,iostat=iostatus)
    end do !i
    read(33,*,iostat=iostatus)
    do i=1,cont2
        cont=cont+1
        read(33,*,iostat=iostatus)teste1,coord(cont,1),coord(cont,2)
    end do !i
    close(33)

    !ELIST
    cont=0
    iostatus=0
    open(44,file='ELIST.lis',status='unknown')
    do while(iostatus==0)
        read(44,*,iostat=iostatus)teste1      
    end do
    close(44)

    read(teste1(1:6),*)cont
    allocate(conec(cont,9))
    conec=0.0
    open(44,file='ELIST.lis',status='unknown')
    read(44,*,iostat=iostatus)
    read(44,*,iostat=iostatus)
    read(44,*,iostat=iostatus)
    cont1=cont/20
    cont2=MOD(cont,20)
    cont=0
    do i=1,cont1
        read(44,*,iostat=iostatus)
        read(44,*,iostat=iostatus)
        do ih=1,20
            cont=cont+1                                 
            read(44,*,iostat=iostatus)teste1,conec(cont,7),teste1,teste1,teste1,teste1,conec(cont,1),&
            conec(cont,3),conec(cont,6),teste1,conec(cont,2),conec(cont,5),teste1,conec(cont,4)
        end do !ih
        read(44,*,iostat=iostatus)
    end do !i
    read(44,*,iostat=iostatus)
    read(44,*,iostat=iostatus)
    do i=1,cont2
        cont=cont+1
        read(44,*,iostat=iostatus)teste1,conec(cont,7),teste1,teste1,teste1,teste1,conec(cont,1),&
        conec(cont,3),conec(cont,6),teste1,conec(cont,2),conec(cont,5),teste1,conec(cont,4)
    end do !i
    close(44)

    !DLIST
    cont=0
    iostatus=0
    open(55,file='DLIST.lis',status='unknown')
    read(55,*,iostat=iostatus)
    read(55,*,iostat=iostatus)
    read(55,*,iostat=iostatus)
    read(55,*,iostat=iostatus)
    do while(iostatus==0)
        read(55,*,iostat=iostatus)teste1
        if(iostatus==0) then
            if(teste1/='NODE  ') cont=cont+1
        end if
    end do
    close(55)

    allocate(dirvinc(cont,3))
    dirvinc=0.0
    open(55,file='DLIST.lis',status='unknown')
    read(55,*,iostat=iostatus)
    read(55,*,iostat=iostatus)
    read(55,*,iostat=iostatus)
    read(55,*,iostat=iostatus)
    cont1=cont/20
    cont2=MOD(cont,20)
    cont=0
    do i=1,cont1
        read(55,*,iostat=iostatus)
        do ih=1,20
            cont=cont+1                                 
            read(55,*,iostat=iostatus)dirvinc(cont,1),teste1,dirvinc(cont,3)
            if(teste1=='UX    ') then
                dirvinc(cont,2)=1.0
            else
                dirvinc(cont,2)=2.0
            end if
        end do !ih
        read(55,*,iostat=iostatus)
    end do !i 
    read(55,*,iostat=iostatus)
    do i=1,cont2
        cont=cont+1                                 
            read(55,*,iostat=iostatus)dirvinc(cont,1),teste1,dirvinc(cont,3)
            if(teste1=='UX    ') then
                dirvinc(cont,2)=1.0
            else
                dirvinc(cont,2)=2.0
            end if
    end do !i
    close(55)

    !FLIST
    cont=0
    iostatus=0
    open(66,file='FLIST.lis',status='unknown')
    read(66,*,iostat=iostatus)
    read(66,*,iostat=iostatus)
    read(66,*,iostat=iostatus)
    read(66,*,iostat=iostatus)
    do while(iostatus==0)
        read(66,*,iostat=iostatus)teste1
        if(iostatus==0) then
            if(teste1/='NODE  ') cont=cont+1
        end if
    end do
    close(66)

    allocate(dircar(cont,3))
    dircar=0.0
    open(66,file='FLIST.lis',status='unknown')
    read(66,*,iostat=iostatus)
    read(66,*,iostat=iostatus)
    read(66,*,iostat=iostatus)
    read(66,*,iostat=iostatus)
    cont1=cont/20
    cont2=MOD(cont,20)
    cont=0
    do i=1,cont1
        read(66,*,iostat=iostatus)
        do ih=1,20
            cont=cont+1                                 
            read(66,*,iostat=iostatus)dircar(cont,1),teste1,dircar(cont,3)
            if(teste1=='FX    ') then
                dircar(cont,2)=1.0
            else
                dircar(cont,2)=2.0
            end if
        end do !ih
        read(66,*,iostat=iostatus)
    end do !i 
    read(66,*,iostat=iostatus)
    do i=1,cont2
        cont=cont+1                                 
            read(66,*,iostat=iostatus)dircar(cont,1),teste1,dircar(cont,3)
            if(teste1=='FX    ') then
                dircar(cont,2)=1.0
            else
                dircar(cont,2)=2.0
            end if
    end do !i
    close(66)
    
    !SFELI
    cont=0
    iostatus=0
    open(77,file='SFELI.lis',status='unknown')
    do while(iostatus==0)
        if(iostatus==0) then
            read(77,*,iostat=iostatus)teste1
            if(teste1=='ELEMEN') then
                do i=1,10
                    read(77,*,iostat=iostatus)cont1,cont2,teste2,value
                    if(cont2==1) then
                        conec(cont1,8)=value
                    else
                        conec(cont1,9)=value
                    end if
                    read(77,*,iostat=iostatus)
                end do !i
            end if
        end if
    end do
    close(77)

    !DADOS_FIBRAS
    cont=0
    iostatus=0
    open(88,file='dados_fibras.txt',status='unknown')
    read(88,*,iostat=iostatus)
    read(88,*,iostat=iostatus)cont
    allocate(nosfibras(cont,2))
    read(88,*,iostat=iostatus)
    do i=1,cont
        read(88,*,iostat=iostatus)nosfibras(i,1),nosfibras(i,2)
    end do
    read(88,*,iostat=iostatus)
    read(88,*,iostat=iostatus)cont
    allocate(conecfibras(cont,2))
    read(88,*,iostat=iostatus)
    do i=1,cont
        read(88,*,iostat=iostatus)conecfibras(i,1),conecfibras(i,2)
    end do
    read(88,*,iostat=iostatus)
    allocate(propfibras(2))
    read(88,*,iostat=iostatus)propfibras(1),propfibras(2)

    close(88)

    end=secnds(start)

    write(*,'(x,a)',advance='no')'LEITURA DE DADOS CONCLUIDA EM '
    write(*,'(f9.5)',advance='no')end
    write(*,'(a)')' SEGUNDOS.'

    write(*,*)' '
    
end subroutine InputData