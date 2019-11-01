subroutine Cplusplus()
    use variables

    implicit none

    integer::i,cont
    real(8)::writestart, writeend

    call cpu_time(writestart)

    write(*,*)'Iniciando exportação de dados para leitura no C++.'

    open(101,file='cplusplus.txt',status='unknown')
    
    write(101,'(a)')'#Materials'
    write(101,*)size(materials,dim=1)
    write(101,'(a)')'Index Young Poisson Density'
    cont=-1
    do i=1,size(materials,dim=1)
        cont = cont+1
        write(101,*)cont, materials(i,1), materials(i,2), 0.0
    end do !i
    
    write(101,'(a)')'#Nodes'
    write(101,*)size(coord,dim=1)
    write(101,'(a)')'Index Coordinates(x1, x2)'
    cont=-1
    do i=1,size(coord,dim=1)
        cont = cont+1
        write(101,*)cont, coord(i,1), coord(i,2)
    end do !i

    write(101,'(a)')'#Elements'
    write(101,*)size(conec,dim=1)
    write(101,'(a)')'Index NodesConection(6) Material Thickness ShapeForce(2)'
    cont = -1
    do i=1,size(conec,dim=1)
        cont = cont+1
        write(101,*)cont, int(conec(i,3)-1), int(conec(i,6)-1), int(conec(i,1)-1), int(conec(i,5)-1),&
        int(conec(i,4)-1), int(conec(i,2)-1), int(conec(i,7)-1), (materials(int(conec(i,7)),3)), conec(i,8), conec(i,9)
    end do !i

    write(101,'(a)')'#NeumannConditions'
    write(101,*)size(dircar,dim=1)
    write(101,'(a)')'Node Direction Value'
    cont = -1
    do i=1,size(dircar,dim=1)
        cont = cont+1
        write(101,*)int(dircar(i,1)-1), int(dircar(i,2)-1), (dircar(i,3))
    end do !i

    write(101,'(a)')'#DirichletConditions'
    write(101,*)size(dirvinc,dim=1)
    write(101,'(a)')'Node Direction Value'
    cont = -1
    do i=1,size(dirvinc,dim=1)
        cont = cont+1
        write(101,*)int(dirvinc(i,1)-1), int(dirvinc(i,2)-1), (dirvinc(i,3))
    end do !i

    call cpu_time(writeend)
    write(*,*)'Exportação do arquivo .txt realizada em', (writeend-writestart),' segundos.'
    write(*,*)' '

end subroutine Cplusplus
