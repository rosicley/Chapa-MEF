subroutine OutputData()
    use variables

    implicit none

    integer::i,nnodes,nelems,nnodesfibras,nelemfibras
    real(4)::start, end


    start=secnds(0.0)

    write(*,'(x,a)')'INICIANDO EXPORTACAO DO ARQUIVO PARA O ACADVIEW.'

    nnodes=size(coord,dim=1)
    nelems=size(conec,dim=1)
    nnodesfibras=size(nosfibras,dim=1)
    nelemfibras=size(conecfibras,dim=1)

    open(99,file='saida.ogl',status='unknown')
    
    write(99,'(a)')'#'
    write(99,100)nnodes+nnodesfibras,nelems+nelemfibras,6
    
    write(99,'(a)')'#'
    do i=1,nnodes
        write(99,200)coord(i,1),coord(i,2),0,0,0,0
    end do !i
    do i=1,nelemfibras
        write(99,200)nosfibras(conecfibras(i,1),1),nosfibras(conecfibras(i,1),2),0,0,0,0
        write(99,200)nosfibras(conecfibras(i,2),1),nosfibras(conecfibras(i,2),2),0,0,0,0
    end do !i

    write(99,'(a)')'#'
    do i=1,nelems
        write(99,300)2,2,int(conec(i,1)),int(conec(i,2)),int(conec(i,3)),&
        int(conec(i,4)),int(conec(i,5)),int(conec(i,6)),0
    end do !i
    do i=1,nelemfibras
        write(99,301)1,1,nnodes+2*(i-1)+1,nnodes+2*(i-1)+2,1
    end do !i

    write(99,'(a)')'#'
    write(99,'(a)')'desloc. x'
    do i=1,nnodes
        write(99,400)sol(2*(i-1)+1),sol(2*(i-1)+2),0,sol(2*(i-1)+1)
    end do !i
    do i=1,nelemfibras
        write(99,400)deslfibras(conecfibras(i,1),1),deslfibras(conecfibras(i,1),2),0,deslfibras(conecfibras(i,1),1)
        write(99,400)deslfibras(conecfibras(i,2),1),deslfibras(conecfibras(i,2),2),0,deslfibras(conecfibras(i,2),1) 
    end do !i
   
    write(99,'(a)')'#'
    write(99,'(a)')'desloc. y'
    do i=1,nnodes
        write(99,400)sol(2*(i-1)+1),sol(2*(i-1)+2),0,sol(2*(i-1)+2)
    end do !i
    do i=1,nelemfibras
        write(99,400)deslfibras(conecfibras(i,1),1),deslfibras(conecfibras(i,1),2),0,deslfibras(conecfibras(i,1),2)
        write(99,400)deslfibras(conecfibras(i,2),1),deslfibras(conecfibras(i,2),2),0,deslfibras(conecfibras(i,2),2) 
    end do !i

    write(99,'(a)')'#'
    write(99,'(a)')'Sigma X'
    do i=1,nnodes
        write(99,400)sol(2*(i-1)+1),sol(2*(i-1)+2),0,nodestress(i,1)
    end do !i
    do i=1,nelemfibras
        write(99,401)deslfibras(conecfibras(i,1),1),deslfibras(conecfibras(i,1),2),0,0
        write(99,401)deslfibras(conecfibras(i,2),1),deslfibras(conecfibras(i,2),2),0,0
    end do !i

    write(99,'(a)')'#'
    write(99,'(a)')'Sigma Y'
    do i=1,nnodes
        write(99,400)sol(2*(i-1)+1),sol(2*(i-1)+2),0,nodestress(i,2)
    end do !i
    do i=1,nelemfibras
        write(99,401)deslfibras(conecfibras(i,1),1),deslfibras(conecfibras(i,1),2),0,0
        write(99,401)deslfibras(conecfibras(i,2),1),deslfibras(conecfibras(i,2),2),0,0
    end do !i

    write(99,'(a)')'#'
    write(99,'(a)')'Tal XY'
    do i=1,nnodes
        write(99,400)sol(2*(i-1)+1),sol(2*(i-1)+2),0,nodestress(i,3)
    end do !i
    do i=1,nelemfibras
        write(99,401)deslfibras(conecfibras(i,1),1),deslfibras(conecfibras(i,1),2),0,0
        write(99,401)deslfibras(conecfibras(i,2),1),deslfibras(conecfibras(i,2),2),0,0
    end do !i

    write(99,'(a)')'#'
    write(99,'(a)')'Normal Fibras'
    do i=1,nnodes
        write(99,401)sol(2*(i-1)+1),sol(2*(i-1)+2),0,0
    end do !i
    do i=1,nelemfibras
        write(99,400)deslfibras(conecfibras(i,1),1),deslfibras(conecfibras(i,1),2),0,normalfibras(i)
        write(99,400)deslfibras(conecfibras(i,2),1),deslfibras(conecfibras(i,2),2),0,normalfibras(i)
    end do !i

    end=secnds(start)

    close(99)

    write(*,'(x,a)',advance='no')'ARQUIVO .OGL EXPORTADO EM '
    write(*,'(f9.5)',advance='no')end
    write(*,'(a)')' SEGUNDOS.'
    write(*,*)' '

    100 format(2x,i10,x,i10,x,i10)
    200 format(x,E16.9,x,E16.9,x,i1,x,i1,x,i1,x,i1)
    !200 format(x,E16.9,x,E16.9,x,f3.1,x,f3.1,x,f3.1,x,f3.1)
    300 format(x,i1,x,i1,x,i10,x,i10,x,i10,x,i10,x,i10,x,i10,x,i1)
    301 format(x,i1,x,i1,x,i10,x,i10,x,i1)
    400 format(x,E16.6,x,E16.6,x,i1,x,E16.6)
    401 format(x,E16.6,x,E16.6,x,i1,x,i1)
    !400 format(x,E16.6,x,E16.6,x,f3.1,x,E16.6)
    !401 format(x,E16.6,x,E16.6,x,f3.1,x,f3.1)

end subroutine OutputData

