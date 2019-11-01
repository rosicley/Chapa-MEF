module variables
    implicit none

    real(8),allocatable::coord(:,:),materials(:,:),conec(:,:),dirvinc(:,:),dircar(:,:),fglobal(:),nosfibras(:,:), &
                         propfibras(:),sol(:),nodestress(:,:),conecnos(:,:),deslfibras(:,:),normalfibras(:)
    character(3)::ep
    integer,allocatable::conecfibras(:,:)
    integer::numberOfThreads






end module variables

