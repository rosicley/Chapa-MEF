subroutine Kernel()
    use sparse
    use variables
    USE omp_lib
    implicit none
   
    integer::i,iel,ih,in,ino,id,j,k,info,nnodes,nelems,ndircar,ndirvinc,no,dir,nnodesfibras,nelemfibras,cont,&
             elem,chapaa(2),auxv
    real(8)::ksi,eta,young,poisson,esp,qx,qy,D11,D22,D12,D33,dxdksi,dxdeta,dydksi,dydeta,peso,det,dksidx,dksidy, &
             detadx,detady,val,dudksi,dudeta,dvdksi,dvdeta,epsx,epsy,gamaxy,sigmx,sigmy,talxy,xc,yc,raio1,raio2, &
             raio3,areaf,youngf,xi,yi,xf,yf,comp,cos,sen
    real(4)::start, end
    real(8)::fi(6),auxMFI(6), auxMDFIDKSI(6),auxMDFIDETA(6)
    real(8)::ksieta(6,2),mat(6,6),coef(6,6),hammer(3,12),MFI(6,12),MdFIdKSI(6,12),MdFIdETA(6,12)
    real(8),allocatable::klocal(:,:),flocal(:),MXX(:,:),MXY(:,:),MYX(:,:),MYY(:,:),cx(:),cy(:),jac(:,:),ijac(:,:), &
                         MDX(:,:),MDY(:,:),MD(:,:),aux1(:,:),aux2(:,:),aux3(:,:),aux4(:,:),MdFIdKSIP(:,:), &
                         MdFIdETAP(:,:),uel(:,:),dataelems(:,:),normf(:),nos(:,:),res(:),kflocal(:,:),&
                         kfglobal(:,:),rotac(:,:),mfibra(:,:),fforma(:),kfc(:,:),dxy(:,:),&
                         uglobal(:),ulocal(:),Faux(:)
    integer,allocatable::ipiv(:),indexes(:),chapa(:)
    integer, dimension(12)::indexes_
    type(sparse_matrix)::spmat

    start = secnds(0.0)
    write(*,'(x,a)')'INICIANDO DETERMINACAO DAS FUNCOES DE FORMAS E DERIVADAS.'
    
    ksieta(1,1) = 0.0
    ksieta(1,2) = 0.0
    ksieta(2,1) = 0.5
    ksieta(2,2) = 0.0
    ksieta(3,1) = 1.0
    ksieta(3,2) = 0.0
    ksieta(4,1) = 0.0
    ksieta(4,2) = 0.5
    ksieta(5,1) = 0.5
    ksieta(5,2) = 0.5
    ksieta(6,1) = 0.0
    ksieta(6,2) = 1.0

    do i=1,6
        ksi = ksieta(i,1)
        eta = ksieta(i,2)
        fi = FIfunction(ksi,eta)
        
        do j=1,6
            mat(j,i)=fi(j)         
        end do !j
    end do !i

    allocate(ipiv(6))

    coef = identityf(6)  
    call DGESV(6,6,mat,6,ipiv,coef,6,info) !resolvendo o sistema

    hammer(1,1)=0.501426509658179
    hammer(1,2)=0.249286745170910
    hammer(1,3)=0.249286745170910
    hammer(1,4)=0.873821971016996
    hammer(1,5)=0.063089014491502
    hammer(1,6)=0.063089014491502
    hammer(1,7)=0.053145049844816
    hammer(1,8)=0.310352451033785
    hammer(1,9)=0.636502499121399
    hammer(1,10)=0.310352451033785
    hammer(1,11)=0.636502499121399
    hammer(1,12)=0.053145049844816

    hammer(2,1)=0.249286745170910
    hammer(2,2)=0.249286745170910
    hammer(2,3)=0.501426509658179
    hammer(2,4)=0.063089014491502
    hammer(2,5)=0.063089014491502
    hammer(2,6)=0.873821971016996
    hammer(2,7)=0.310352451033785
    hammer(2,8)=0.636502499121399
    hammer(2,9)=0.053145049844816
    hammer(2,10)=0.053145049844816
    hammer(2,11)=0.310352451033785
    hammer(2,12)=0.636502499121399

    hammer(3,1)=0.05839313786319
    hammer(3,2)=0.05839313786319
    hammer(3,3)=0.05839313786319
    hammer(3,4)=0.025422453185104
    hammer(3,5)=0.025422453185104
    hammer(3,6)=0.025422453185104
    hammer(3,7)=0.041425537809187
    hammer(3,8)=0.041425537809187
    hammer(3,9)=0.041425537809187
    hammer(3,10)=0.041425537809187
    hammer(3,11)=0.041425537809187
    hammer(3,12)=0.041425537809187

    do i=1,12
        ksi=hammer(1,i)
        eta=hammer(2,i)
        auxMFI=matmul(coef,FIfunction(ksi,eta))
        auxMDFIDKSI=matmul(coef,dFIdKSIfunction(ksi,eta))
        auxMDFIDETA=matmul(coef,dFIdETAfunction(ksi,eta))

        do j=1,6
            MFI(j,i)=auxMFI(j)
            MdFIdKSI(j,i)=auxMDFIDKSI(j)
            MdFIdETA(j,i)=auxMDFIDETA(j)
        end do !j
    end do !i

    allocate(MdFIdKSIP(6,6),MdFIdETAP(6,6))
    
    do i=1,6
        ksi=ksieta(i,1)
        eta=ksieta(i,2)
        auxMDFIDKSI=matmul(coef,dFIdKSIfunction(ksi,eta))
        auxMDFIDETA=matmul(coef,dFIdETAfunction(ksi,eta))

        do j=1,6
            MdFIdKSIP(j,i)=auxMDFIDKSI(j)
            MdFIdETAP(j,i)=auxMDFIDETA(j)

        end do !j
    end do !i

    !!dirvinc e dircar(nó, direção (1=x e 2=y), valor)
    
    nnodes=size(coord,dim=1)
    nelems=size(conec,dim=1)
    ndircar=size(dircar,dim=1)
    ndirvinc=size(dirvinc,dim=1)
    nnodesfibras=size(nosfibras,dim=1)
    nelemfibras=size(conecfibras,dim=1)

    if(nnodesfibras>1) then
        call prepare_to_use(spmat, (2*nnodes), (nelems*78+nelemfibras*576))
    else
        call prepare_to_use(spmat, (2*nnodes), (nelems*78))
    end if

    end = secnds(start)

    write(*,'(x,a)',advance='no')'OPERACAO REALIZADA EM '
    write(*,'(f9.5)',advance='no')end
    write(*,'(a)')' SEGUNDOS.'

    write(*,*)' '
    
    CALL omp_set_num_threads(numberOfThreads)

    
    if(nnodesfibras>1) then
        start = secnds(0.0)
        write(*,'(x,a)')'INICIANDO CONTRIBUICAO DAS FIBRAS NO DOMINIO DO SOLIDO.'


        allocate(dataelems(nelems,3),normf(2),nos(6,2),res(3),conecnos(nnodesfibras,3))
        conecnos=0.0d0
        nos=0.0d0
        res=0.0d0

        youngf=propfibras(1)
        areaf=propfibras(2)

        allocate(chapa(2),kflocal(4,4),kfglobal(4,4),rotac(4,4),mfibra(24,4),fforma(6),&
                kfc(24,24),indexes(24))

        do i=1,nelems
            xc=(coord(int(conec(i,1)),1)+coord(int(conec(i,3)),1)+coord(int(conec(i,6)),1))/3
            yc=(coord(int(conec(i,1)),2)+coord(int(conec(i,3)),2)+coord(int(conec(i,6)),2))/3

            raio1=norm2((/xc-coord(int(conec(i,1)),1),yc-coord(int(conec(i,1)),2)/))
            raio2=norm2((/xc-coord(int(conec(i,3)),1),yc-coord(int(conec(i,3)),2)/))
            raio3=norm2((/xc-coord(int(conec(i,6)),1),yc-coord(int(conec(i,6)),2)/))

            dataelems(i,:)=(/xc,yc,max(raio1,raio2,raio2)/)
        end do !i

        !$OMP PARALLEL PRIVATE(normf,nos,res)
        !$OMP DO
        do i=1,nnodesfibras
            do iel=1,nelems
                normf(1)=nosfibras(i,1)-dataelems(iel,1)
                normf(2)=nosfibras(i,2)-dataelems(iel,2)
                

                if(norm2(normf)<=2*dataelems(iel,3)) then
                    do j=1,6
                        !nos(j,:)=(/coord(int(conec(iel,j)),1),coord(int(conec(iel,j)),2)/)
                        nos(j,:)=coord(int(conec(iel,j)),:)
                    end do !j
                    
                    res=fksieta(nos,nosfibras(i,:))

                    if((res(1)>=0.0d0).and.(res(1)<=1.0d0).and.&
                       (res(2)>=0.0d0).and.(res(2)<=1.0d0).and.&
                       (res(3)>=0.0d0).and.(res(3)<=1.0d0)) then

                        conecnos(i,:)=(/iel*1.0d0,res(1),res(2)/)

                        EXIT

                    end if
                end if
            end do !iel
        end do !i
        !$OMP END DO
        !$OMP END PARALLEL
    
          
        
        !$OMP PARALLEL PRIVATE(chapa,xi,yi,xf,yf,comp,cos,sen,kflocal,rotac,kfglobal,mfibra,ksi,eta,fforma,kfc,indexes,cont)
        !$OMP DO
        do iel=1,nelemfibras
            chapa=0
            chapa(1)=int(conecnos(conecfibras(iel,1),1))
            chapa(2)=int(conecnos(conecfibras(iel,2),1))

            if(chapa(1)*chapa(2)>0) then
                xi=nosfibras(conecfibras(iel,1),1)
                yi=nosfibras(conecfibras(iel,1),2)
                xf=nosfibras(conecfibras(iel,2),1)
                yf=nosfibras(conecfibras(iel,2),2)

                comp=((xf-xi)*(xf-xi)+(yf-yi)*(yf-yi))**0.50d0

                cos=(xf-xi)/comp
                sen=(yf-yi)/comp

                kflocal(1,:)=youngf*areaf/comp*(/1.0d0,0.0d0,-1.0d0,0.0d0/)
                kflocal(2,:)=youngf*areaf/comp*(/0.0d0,0.0d0,0.0d0,0.0d0/)
                kflocal(3,:)=youngf*areaf/comp*(/-1.0d0,0.0d0,1.0d0,0.0d0/)
                kflocal(4,:)=youngf*areaf/comp*(/0.0d0,0.0d0,0.0d0,0.0d0/)

                rotac(1,:)=(/cos,-sen,0.0d0, 0.0d0/)
                rotac(2,:)=(/sen,cos,0.0d0, 0.0d0/)
                rotac(3,:)=(/0.0d0, 0.0d0,cos,-sen/)
                rotac(4,:)=(/0.0d0, 0.0d0,sen,cos/)

                kfglobal=matmul(matmul(rotac,kflocal),transpose(rotac))

                mfibra=0.0d0
                ksi=conecnos(conecfibras(iel,1),2)
                eta=conecnos(conecfibras(iel,1),3)
                fforma=matmul(coef,FIfunction(ksi,eta))

                do i=1,6
                    mfibra(2*(i-1)+1,1)=fforma(i)
                    mfibra(2*(i-1)+2,2)=fforma(i)
                end do !i

                ksi=conecnos(conecfibras(iel,2),2)
                eta=conecnos(conecfibras(iel,2),3)
                fforma=matmul(coef,FIfunction(ksi,eta))

                do i=1,6
                    mfibra(12+2*(i-1)+1,3)=fforma(i)
                    mfibra(12+2*(i-1)+2,4)=fforma(i)
                end do !i

                kfc=matmul(matmul(mfibra,kfglobal),transpose(mfibra))
                !vindexes=0.0d0
                indexes=0

                cont=0
                do i=1,2
                    do j=1,6
                        do k=1,2
                            cont=cont+1
                            indexes(cont)=2*(conec(int(chapa(i)),j)-1)+k
                        end do !k
                    end do !j
                end do !i

                !$OMP CRITICAL
                call ADD_FULL_MATRIX(spmat,kfc,indexes,24)
                !$OMP END CRITICAL

            end if
        end do !iel 
        !$OMP END DO 
        !$OMP END PARALLEL


        end=secnds(start)

        write(*,'(x,a)',advance='no')'CONTRIBUICAO DAS FIBRAS REALIZADA EM '
        write(*,'(f9.5)',advance='no')end
        write(*,'(a)')' SEGUNDOS.'
        write(*,*)' '
    end if

    start=secnds(0.0)
    write(*,'(x,a)')'INICIANDO CONTRIBUICAO DOS ELEMENTOS DE CHAPA E MONTAGEM DO VETOR DE FORCAS NODAIS NO DOMINIO DO SOLIDO.'

    
    allocate(MXX(2,2),MXY(2,2),MYX(2,2),MYY(2,2),cx(6),cy(6),jac(2,2),ijac(2,2),MDX(4,2),MDY(4,2), &
             MD(12,4),klocal(12,12),flocal(12),aux1(4,4),aux2(4,4),aux3(4,4),aux4(4,4),fglobal(2*nnodes)) !indexes_ ja foi alocado no começo, ver para paralelizar

    MD=0.0d0
    fglobal=0.0d0
    

    !Montagem da Matriz de Rigidez e Vetor de Forças Nodais Equivalentes
    !$OMP PARALLEL FIRSTPRIVATE(flocal,klocal,young,poisson,esp,qx,qy,D11,D22,D12,D33,MXX,MXY,MYX,MYY,cx,cy,peso,aux2,aux3,aux4,&
    !$OMP dxdksi,dxdeta,dydksi,dydeta,jac,det,ijac,dksidx,dksidy,detadx,detady,MDX,MDY,MD,aux1,cont,indexes_)
    !$OMP DO
    do iel=1,nelems
        flocal=0.0
        klocal=0.0
        young=materials(int(conec(iel,7)),1)
        poisson=materials(int(conec(iel,7)),2)
        esp=materials(int(conec(iel,7)),3)
        qx=conec(iel,8)
        qy=conec(iel,9)
        if(ep=='EPT') then
            D11=young/(1-poisson*poisson)
            D22=young/(1-poisson*poisson)
            D12=poisson*young/(1-poisson*poisson)
        else 
            D11=(1-poisson)*young/(1+poisson)/(1-2*poisson)
            D22=(1-poisson)*young/(1+poisson)/(1-2*poisson)
            D12=poisson*young/(1+poisson)/(1-2*poisson)
        end if

        D33=young/(1+poisson)
                
        MXX(1,:)=(/D11,0.0d0/)
        MXX(2,:)=(/0.0d0,0.5d0*D33/)
       
        MXY(1,:)=(/0.0d0,D12/)
        MXY(2,:)=(/0.5d0*D33,0.0d0/)

        MYX(1,:)=(/0.0d0,0.5d0*D33/)
        MYX(2,:)=(/D12,0.0d0/)

        MYY(1,:)=(/0.5d0*D33,0.0d0/)
        MYY(2,:)=(/0.0d0,D22/)

        do i=1,6
            cx(i)=coord(int(conec(iel,i)),1)
            cy(i)=coord(int(conec(iel,i)),2)
        end do !i
        
        do ih=1,12 !12 é o número de hammer
            peso=hammer(3,ih)
            dxdksi=0.0
            dxdeta=0.0
            dydksi=0.0
            dydeta=0.0

            do j=1,6
                dxdksi=dxdksi+cx(j)*MdFIdKSI(j,ih)
                dxdeta=dxdeta+cx(j)*MdFIdETA(j,ih)
                dydksi=dydksi+cy(j)*MdFIdKSI(j,ih)
                dydeta=dydeta+cy(j)*MdFIdETA(j,ih)
            end do !j

            jac(1,:)=(/dxdksi,dxdeta/)
            jac(2,:)=(/dydksi,dydeta/)
            det=dxdksi*dydeta-dxdeta*dydksi

            ijac=matinv2(jac) 
            
            dksidx=ijac(1,1)
            dksidy=ijac(1,2)
            detadx=ijac(2,1)
            detady=ijac(2,2)

            MDX(1,:)=(/dksidx,0.0d0/)
            MDX(2,:)=(/0.0d0,dksidx/)
            MDX(3,:)=(/detadx,0.0d0/)
            MDX(4,:)=(/0.0d0,detadx/)

            MDY(1,:)=(/dksidy,0.0d0/)
            MDY(2,:)=(/0.0d0,dksidy/)
            MDY(3,:)=(/detady,0.0d0/)
            MDY(4,:)=(/0.0d0,detady/)

            do i=1,6
                MD(2*(i-1)+1,1)=MdFIdKSI(i,ih)
                MD(2*(i-1)+1,3)=MdFIdETA(i,ih)
                MD(2*(i-1)+2,2)=MdFIdKSI(i,ih)
                MD(2*(i-1)+2,4)=MdFIdETA(i,ih)
            end do !i

            do i=1,6
                flocal(2*(i-1)+1)=flocal(2*(i-1)+1)+qx*MFI(i,ih)*det*peso
                flocal(2*(i-1)+2)=flocal(2*(i-1)+2)+qy*MFI(i,ih)*det*peso
            end do !i

            aux1=matmul(matmul(MDX,MXX),transpose(MDX))
            aux2=matmul(matmul(MDX,MXY),transpose(MDY))
            aux3=matmul(matmul(MDY,MYX),transpose(MDX))
            aux4=matmul(matmul(MDY,MYY),transpose(MDY))

            klocal=klocal+matmul(matmul(MD,(aux1+aux2+aux3+aux4)),transpose(MD))*esp*det*peso      
        end do !ih

        cont = 0

        do in=1,6
            do id=1,2
                cont = cont +1
                indexes_(cont)=2*(int(conec(iel,in))-1)+id
                !$OMP CRITICAL
                fglobal(2*(int(conec(iel,in))-1)+id)=fglobal(2*(int(conec(iel,in))-1)+id)+flocal(2*(in-1)+id)            
                !$OMP END CRITICAL
            end do !id
        end do !in

        !$OMP CRITICAL
        call add_matrix(spmat, klocal, indexes_, 12)
        !$OMP END CRITICAL

    end do !iel
    !$OMP END DO 
    !$OMP END PARALLEL

    end=secnds(start)

    write(*,'(x,a)',advance='no')'CONTRIBUICAO DOS ELEMENTOS DE CHAPA REALIZADA EM '
    write(*,'(f9.5)',advance='no')end
    write(*,'(a)')' SEGUNDOS.'
    write(*,*)' '


    start=secnds(0.0)
    write(*,'(x,a)')'INICIANDO MONTAGEM DA MATRIZ SPARSA E APLICACAO DAS CONDICOES DE CONTORNO.'


    call assemble_sparse_matrix(spmat, timeit=.false.)
    call optimize_storage(spmat)

    
    !Aplicação das condições de contorno em forças
    do i=1,ndircar
        no=int(dircar(i,1))
        dir=int(dircar(i,2))
        val=dircar(i,3)
        fglobal(2*(no-1)+dir)=fglobal(2*(no-1)+dir)+val
    end do !i

    allocate(Faux(2*nnodes))
    !Aplicação das condições de contorno em deslocamentos
    do i=1,ndirvinc
        auxv=2*(int(dirvinc(i,1))-1)+int(dirvinc(i,2))
        val=dirvinc(i,3)
        call GET_ROW2(spmat,auxv,Faux)
        fglobal=fglobal-val*Faux

        CALL SET_VALUE_TO_ROW(spmat,auxv,0.d0)
        CALL SET_VALUE_TO_COL(spmat,auxv,0.d0)
        CALL SET_VALUE_IN_TERM(spmat,auxv,auxv,1.d0)

        fglobal(auxv)=val
    end do !i

    end=secnds(start)
    write(*,'(x,a)',advance='no')'OPERACAO REALIZADA EM '
    write(*,'(f9.5)',advance='no')end
    write(*,'(a)')' SEGUNDOS.'
    write(*,*)' '

    write(*,'(x,a)')'INICIANDO SOLUCAO DO SISTEMA LINEAR.'
    start = secnds(0.0)

    deallocate(ipiv)
    allocate(sol(2*nnodes))


    sol=fglobal
    !call DPOSV('U',2*nnodes,1,kglobal,2*nnodes,sol,2*nnodes,info)
    call solve_system_of_equation(spmat, sol)
         
    end=secnds(start)
    write(*,'(x,a)',advance='no')'SISTEMA RESOLVIDO EM '
    write(*,'(f9.5)',advance='no')end
    write(*,'(a)')' SEGUNDOS.'
    write(*,*)' '

    start=secnds(0.0)
    write(*,'(x,a)')'INICIANDO CALCUDO DAS TENSOES NODAIS E NORMAIS NAS FIBRAS.'

    !Tensões
    allocate(nodestress(nnodes,4),uel(6,2))
    nodestress=0.d0

    do iel=1,nelems
        uel=0.d0
        do ino=1,6
            uel(ino,1)=sol(2*(int(conec(iel,ino))-1)+1)
            uel(ino,2)=sol(2*(int(conec(iel,ino))-1)+2)
        end do !ino
        young=materials(int(conec(iel,7)),1)
        poisson=materials(int(conec(iel,7)),2)
        esp=materials(int(conec(iel,7)),3)

        if(ep=='EPT') then
            D11=young/(1-poisson*poisson)
            D22=young/(1-poisson*poisson)
            D12=poisson*young/(1-poisson*poisson)
        else
            D11=(1-poisson)*young/(1+poisson)/(1-2*poisson)
            D22=(1-poisson)*young/(1+poisson)/(1-2*poisson)
            D12=poisson*young/(1+poisson)/(1-2*poisson)
        end if
        D33=young/(1+poisson)

        do i=1,6
            cx(i)=coord(int(conec(iel,i)),1)
            cy(i)=coord(int(conec(iel,i)),2)
        end do !i
        
        do ino=1,6
            dxdksi=0.d0
            dxdeta=0.d0
            dydksi=0.d0
            dydeta=0.d0
            do j=1,6
                dxdksi=dxdksi+cx(j)*MdFIdKSIP(j,ino)
                dxdeta=dxdeta+cx(j)*MdFIdETAP(j,ino)
                dydksi=dydksi+cy(j)*MdFIdKSIP(j,ino)
                dydeta=dydeta+cy(j)*MdFIdETAP(j,ino)
            end do !j
            dudksi=0.d0
            dudeta=0.d0
            dvdksi=0.d0
            dvdeta=0.d0

            do i=1,6
                dudksi=dudksi+uel(i,1)*MdFIdKSIP(i,ino)
                dudeta=dudeta+uel(i,1)*MdFIdETAP(i,ino)
                dvdksi=dvdksi+uel(i,2)*MdFIdKSIP(i,ino)
                dvdeta=dvdeta+uel(i,2)*MdFIdETAP(i,ino)
            end do !i

            jac(1,:)=(/dxdksi,dxdeta/)
            jac(2,:)=(/dydksi,dydeta/)
            det=dxdksi*dydeta-dxdeta*dydksi
            
            ijac=matinv2(jac)

            dksidx=ijac(1,1)
            dksidy=ijac(1,2)
            detadx=ijac(2,1)
            detady=ijac(2,2)
            epsx = dudksi*dksidx+dudeta*detadx
            epsy = dvdksi*dksidy+dvdeta*detady
            gamaxy = 0.5d0*(dudksi*dksidy+dudeta*detady+dvdksi*dksidx+dvdeta*detadx)
            sigmx = D11*epsx+D12*epsy
            sigmy = D12*epsx+D22*epsy
            talxy = D33*gamaxy
            
            nodestress(int(conec(iel,ino)),1)=nodestress(int(conec(iel,ino)),1)+sigmx
            nodestress(int(conec(iel,ino)),2)=nodestress(int(conec(iel,ino)),2)+sigmy
            nodestress(int(conec(iel,ino)),3)=nodestress(int(conec(iel,ino)),3)+talxy
            nodestress(int(conec(iel,ino)),4)=nodestress(int(conec(iel,ino)),4)+1.d0
        end do !ino
    end do !iel
    do i=1,3
        do j=1,nnodes
        nodestress(j,i)=nodestress(j,i)/nodestress(j,4)
        end do !j
    end do !i

    allocate(deslfibras(nnodesfibras,2),dxy(6,2))
    deslfibras=0.d0

    do i=1,nnodesfibras
        elem=int(conecnos(i,1))
        if(elem>0) then
            ksi=conecnos(i,2)
            eta=conecnos(i,3)
            fforma=matmul(coef,FIfunction(ksi,eta))

            do ih=1,6
                do j=1,2
                    dxy(ih,j)=sol(2*(int(conec(elem,ih))-1)+j)
                end do !j
            end do !ih

            deslfibras(i,1)=dot_product(fforma,dxy(:,1))
            deslfibras(i,2)=dot_product(fforma,dxy(:,2))
        end if
    end do !i

    
    allocate(normalfibras(nelemfibras),uglobal(4),ulocal(4))

    normalfibras=0.0d0
    uglobal=0.0d0
    ulocal=0.0d0
    youngf=propfibras(1)
    areaf=propfibras(2)

    do i=1,nelemfibras
        chapaa=0
        chapaa(1)=int(conecnos(conecfibras(i,1),1))
        chapaa(2)=int(conecnos(conecfibras(i,2),1))

        if(chapaa(1)*chapaa(2)>0) then
            xi=nosfibras(conecfibras(i,1),1)
            yi=nosfibras(conecfibras(i,1),2)
            xf=nosfibras(conecfibras(i,2),1)
            yf=nosfibras(conecfibras(i,2),2)

            comp=((xf-xi)*(xf-xi)+(yf-yi)*(yf-yi))**0.5

            cos=(xf-xi)/comp
            sen=(yf-yi)/comp

            rotac(1,:)=(/cos,-sen,0.d0, 0.d0/)
            rotac(2,:)=(/sen,cos,0.d0, 0.d0/)
            rotac(3,:)=(/0.d0, 0.d0,cos,-sen/)
            rotac(4,:)=(/0.d0, 0.d0,sen,cos/)

            uglobal=(/deslfibras(conecfibras(i,1),1),deslfibras(conecfibras(i,1),2),&
                      deslfibras(conecfibras(i,2),1),deslfibras(conecfibras(i,2),2)/)

            ulocal=matmul(transpose(rotac),uglobal)
            normalfibras(i)=(ulocal(3)-ulocal(1))*youngf*areaf/comp
        end if
    end do !i

    end=secnds(start)
    write(*,'(x,a)',advance='no')'CALCULO REALIZADO EM '
    write(*,'(f9.5)',advance='no')end
    write(*,'(a)')' SEGUNDOS.'
    write(*,*)' '
    
    contains
        function FIfunction(ksi,eta) result (vec)
            real(8)::vec(6)
            real(8)::ksi, eta
            vec = (/1.0d0, ksi, eta, ksi*ksi, eta*eta, ksi*eta/)
        end function FIfunction

        function dFIdKSIfunction(ksi,eta) result (vec)
            real(8)::vec(6)
            real(8)::ksi, eta
            vec = (/0.0d0, 1.0d0, 0.0d0, 2.0d0*ksi, 0.0d0, eta/)
        end function dFIdKSIfunction

        function dFIdETAfunction(ksi,eta) result (vec)
            real(8)::vec(6)
            real(8)::ksi, eta
            vec = (/0.0d0, 0.0d0, 1.0d0, 0.0d0, 2.0d0*eta, ksi/)
        end function dFIdETAfunction     

        function identityf(number) result (matrix)
            integer::number,i,j
            real(8)::matrix(number,number)
            matrix=0.0
            do i=1,number
                do j=1,number
                    if(i==j) matrix(i,j)=1.0
                end do
            end do

        end function identityf

        function matinv2(A) result(B)
            !! Performs a direct calculation of the inverse of a 2×2 matrix.
            real(8),intent(in) :: A(2,2)   !! Matrix
            real(8):: B(2,2)   !! Inverse matrix
            real(8):: detinv
        
            ! Calculate the inverse determinant of the matrix
            detinv = 1/(A(1,1)*A(2,2) - A(1,2)*A(2,1))
        
            ! Calculate the inverse of the matrix
            B(1,1) = +detinv * A(2,2)
            B(2,1) = -detinv * A(2,1)
            B(1,2) = -detinv * A(1,2)
            B(2,2) = +detinv * A(1,1)
      end function

    function fksieta(nodes,ponto) result (sol1)
        real(8)::nodes(6,2),ponto(2),sol1(3),p1(2),p2(2),p3(2),v13(2),v23(2),m1(2,2), &
                 m(2,2),auxsol(2),ksi0,eta0,ponto0(2),dsol(2),erro
        integer::cont

        p3=nodes(1,:)
        p1=nodes(3,:)
        p2=nodes(6,:)

        v13=p1-p3
        v23=p2-p3

        m1(1,:)=v13
        m1(2,:)=v23

        m=matinv2(transpose(m1))
        auxsol=matmul(m,(ponto-p3))

        erro=1.0d0
        cont=0

        do while((erro>=0.000001).and.(cont<=20))
            ksi0=auxsol(1)
            eta0=auxsol(2)

            ponto0(1)=dot_product(matmul(coef,FIfunction(ksi0,eta0)),nodes(:,1))
            ponto0(2)=dot_product(matmul(coef,FIfunction(ksi0,eta0)),nodes(:,2))

            m1(1,1)=dot_product(matmul(coef,dFIdKSIfunction(ksi0,eta0)),nodes(:,1))
            m1(1,2)=dot_product(matmul(coef,dFIdETAfunction(ksi0,eta0)),nodes(:,1))
            m1(2,1)=dot_product(matmul(coef,dFIdKSIfunction(ksi0,eta0)),nodes(:,2))
            m1(2,2)=dot_product(matmul(coef,dFIdETAfunction(ksi0,eta0)),nodes(:,2))

            m=matinv2(m1)
            dsol=matmul(m,(ponto-ponto0))
            auxsol=auxsol+dsol

            erro=norm2(dsol)
            cont=cont+1
        end do

        sol1=(/auxsol(1),auxsol(2),1.0d0-auxsol(1)-auxsol(2)/)
    end function fksieta

end subroutine Kernel