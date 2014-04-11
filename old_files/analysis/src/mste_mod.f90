module Parameters
  character(20) :: fname
  real (kind=8) :: rc
  integer :: nvec
  integer :: iconf=0, nconf=0
  logical :: ipbc
  integer :: fstat
end module Parameters

module Particles
  integer       :: npart, nprot
  real (kind=8) :: boxsz
  real (kind=8), dimension(:,:), allocatable :: x, p, x_shifted
  logical, dimension(:), allocatable :: isos
  real (kind=8) ,parameter :: mass=939.0
end module Particles

module Statistics
  integer :: pasada
end module Statistics

module mod_mst
  use Particles
  use Parameters
  implicit none
  integer     :: nclus0, nclus1 ! Nro de clusters antes y despues de PBC
  integer     :: mass_inf
  ! Flags
  logical     :: isospin

  ! Datos fisicos

  Real,parameter :: ee=2.81777*0.511006 !charge squared for Coulomb
  Real,parameter :: lambda=20.0
  Real,parameter :: invlambda=1.0/lambda

  !   Parametros Pandha Medium
  Real,parameter :: vr=3088.118
  Real,parameter :: va=2666.647
  Real,parameter :: v0=373.118
  Real,parameter :: mur=1.7468
  Real,parameter :: mua=1.6
  Real,parameter :: mu0=1.5
  Real,parameter :: rc0=4.5

  Real,parameter :: ep0_nn=(1.d0/rc0)*v0*exp(-mu0*rc0)
  Real,parameter :: ep0_np=(1.d0/rc0)*(vr*exp(-mur*rc0) - va*exp(-mua*rc0))

  ! Listas sin sentido fisico
  integer, allocatable    :: lista_part_rem(:) ! Lista de particulas que falta clusterizar
  ! (ex ivwork)
  integer, allocatable    :: lista_part_aux(:) ! Lista de particulas que se clusterizan en 
  ! cada paso (ex ivaux)
!!!!!!!!!
  integer, allocatable    :: index_part_clu(:) ! Indice de particulas que se agregan
  ! a un cluster en un paso (ex ivind)

  integer, allocatable    :: lista_part_ord(:) ! Eso. (ex iaapar)
  integer, allocatable    :: bordes_clusters(:)! Lista de indices de 'bordes' de clusters
  ! (ex iaalim). Debe alocarse de 0 a npart
  ! Cosas con sentido fisico
  integer, allocatable    :: mass_dist(:),mdistt(:)      ! Distribucion de masas
Contains

  !______________________________________________________________________________
  !______________________________________________________________________________
  !Subrutina que calcula fragmentos MST sin PBC
  !______________________________________________________________________________
  !______________________________________________________________________________
  subroutine mst()
    Implicit none

    real            :: x1(3), p1(3),dp(3),pp,d(3), dd, rc2
    real            :: e12
    integer         :: ipart, jpart, kpart ! Iteradores de particulas
    integer         :: iaux1, iaux2, isort1, isort2! Iteradores genericos
    !Se usan para reordenar listas
    integer         :: rempart ! Nro. de particulas no clusterizadas aun
    integer         :: iclus   ! Indice de cluster
    integer         :: part_in_clu  ! Nro. de particulas que se asignan
    ! a un cluster en cada paso.
    integer         :: part_added   ! Nro. de particulas que se agregan
    ! al cluster
    integer         :: ip1,ip2 ! Indices de particulas a comparar
    ! Inicializo cosas
    iclus = 0
    lista_part_aux=0
    lista_part_rem=0
    lista_part_ord=0
    bordes_clusters=0
    index_part_clu=0
    mass_dist=0
    rc=rc0
    rc2=rc*rc

    rempart = npart ! Remaining particles: Particulas que falta clusterizar
    ! hasta ahora
    iclus = 0       ! Numero de clusters hasta ahora (e indice de cluster)

    do ipart = 1, npart
       lista_part_rem(ipart) = ipart ! Lista de particulas que falta 
       ! clusterizar ("remanentes")
    enddo
    ! Go.
    do ipart = 1, npart ! Recorro particulas (por primera vez)
       iclus = iclus + 1 ! Indice de cluster. Inauguro cluster nuevo.
       lista_part_aux=0  ! Lista de particulas que van a formar el cluster
       ! iclus-esimo
       lista_part_aux(1) = lista_part_rem(1) ! Empiezo a clusterizar en el orden
       ! que vienen
       rempart=rempart-1 ! Particulas que falta clusterizar hasta ahora
       if(rempart<0) exit ! Si ya revise todas, salgo.
       ! Voy a buscar particulas que 'clustericen' con ipart (la primera de la lista
       ! de remanentes)
       ! Saco a ipart de la lista de remanentes y agrupo las otras particulas al
       ! principio de esa lista.
       do iaux1 = 1,rempart
          lista_part_rem(iaux1)=lista_part_rem(iaux1+1)
       enddo
       ! Taggeo el final de la lista con un cero
       lista_part_rem(rempart+1) = 0
       part_in_clu = 1 ! Numero de particulas en cluster que construyo en
       ! este ciclo
       ! Ahora recorro todas las particulas que ya estaban en el cluster, para ver si
       ! alguna de las remanentes se pega al iclus-esimo debido a la ultima que se
       ! agrego.
       do jpart = 1, npart ! Las particulas ya clusterizadas estan
          ! ordenadas al principio de lista_part_aux
          if (jpart > part_in_clu) exit ! Salgo cuando ya revise todas
          ip1 = lista_part_aux(jpart) ! Indice de particula ya clusterizada
          x1(:)=x(:,ip1) ! Coordenadas de dicha particula
          p1(:)=p(:,ip1) ! Momentos de dicha particula
          ! Ahora recorro las remanentes y me fijo si alguna se agrega al cluster
          ! iclus-esimo por estar cerca de ip1
          part_added = 0 ! Cuento cuantas particulas se agregan al cluster
          ! iclus-esimo por estar cerca de ip1
          do kpart = 1, rempart ! Particulas todavia no clusterizadas.
             ip2 = lista_part_rem(kpart)
             ! Me fijo si clusterizan (MST, MSTE o MSD difieren solo en el criterio)
             d(:) = x1(:)-x(:,ip2)
             where(d>boxsz/2) d=d-boxsz
             where(d<-boxsz/2) d=d+boxsz
             dd = dot_product(d,d)
             ! MST Fisico: Si estan cerca y son de distinto isospin, se pegan
             !          if(isospin) then
             if(dd<rc2 .and. (isos(ip1).neqv.isos(ip2))) then
                dp(:)=p1(:)-p(:,ip2)
                pp = dot_product(dp,dp)
                e12 = epot(dd,isos(ip1),isos(ip2)) + pp/(4.0*mass)
                if(e12<=0) then       ! Criterio para MSTE
                   rempart = rempart-1 ! Queda una menos por clusterizar
                   part_added = part_added + 1 ! Se agrega una particula mas
                   ! al cluster de ip1
                   part_in_clu = part_in_clu + 1 ! Hay una mas clusterizada overall
                   lista_part_aux(part_in_clu) = ip2 ! Agrego
                   ! la nueva particula del cluster (ip2) al lado de la ultima
                   ! que habia.
                   index_part_clu(part_added) = kpart ! Anoto particulas a sacar de
                   ! la lista de remanentes 
                endif ! clusteriza con MSTE
             endif
             
             !          else ! No chequeo isospin
             !          endif ! Chequeo isospin o no
          enddo ! Todas las particulas que no estaban clusterizadas
          
          ! Ahora la parte delicada: Reordenar las listas.
          ! Tengo que sacar las clusterizadas de la lista de remanentes,reordenar
          ! todas las remanentes (correr a la izquierda) para que ocupen las primeras
          ! rempart posiciones, y rellenar con ceros

          ! Si se agrega una unica particula
          if (part_added==1) then
             ! Corro una posicion a la izquierda todas las particulas que estan a la
             ! derecha de esa y listo.
             do iaux1 = index_part_clu(part_added)+1, npart
                lista_part_rem(iaux1-1)=lista_part_rem(iaux1)
             enddo
             ! Pero si se agrega mas de una, la cosa se complica:
          else if (part_added>1) then
             ! Recorro las que se agregaron. Estan indexadas en index_part_clu
             do iaux1 = 1, part_added-1
                ! Identifico la primera particula a la derecha de una particula agregada
                ! al cluster en este paso.
                isort1 = index_part_clu(iaux1)+1
                ! Identifico la primera particula a la izquierda de la SIGUIENTE particula
                ! agregada.
                isort2 = index_part_clu(iaux1+1)-1
                ! Las particulas entre isort1 e isort2 en lista_part_rem siguen siendo
                ! remanentes y hay que correrlas a la izquierda.
                ! Las recorro...
                do iaux2 = isort1,isort2
                   ! y calculo donde tengo que correrlas:
                   !  "Posicion en la que estaban" - "# de particulas clusterizadas en este paso
                   !   que estaban a la izquierda"
                   lista_part_rem(iaux2-iaux1) = lista_part_rem(iaux2)
                enddo
             enddo  !cluster size > 1
             ! Ademas, pueden quedar particulas sin clusterizar a la derecha de la ultima
             ! que fue cluterizada en este paso.
             isort2 = isort2+2 ! Posicion de la primera remanente a la derecha
             ! de la ultima clusterizada.
             ! Las corro a la izquierda part_added posiciones
             do iaux2 = isort2, npart
                lista_part_rem(iaux2-part_added) = lista_part_rem(iaux2)
             enddo
             ! Y relleno con ceros el resto
             lista_part_rem(rempart+1:npart) = 0
          endif ! cluster size query
       enddo ! Agrandar cluster
       ! Si llego a esta linea es porque acabo de terminar de reconocer un fragmento
       ! nuevo de tamaño part_in_clu. Lo agrego a la ditribucion de masas.
       mass_dist(part_in_clu) = mass_dist(part_in_clu) + 1
       ! Compagino las listas definitivas
       ! 1) Lista ordenada de particulas. Concateno las part_in_clu despues del ultimo
       ! cluster ya identificado:
       iaux1 = bordes_clusters(iclus-1)+1 ! Posicion de la primera particula
       ! del cluster nuevo en la lista
       ! definitiva (lista_part_ord)
       iaux2 = iaux1+part_in_clu-1        ! Posicion de la ultima particula
       ! del cluster nuevo
       lista_part_ord(iaux1:iaux2) = lista_part_aux(1:part_in_clu)
       ! Lista de 'bordes'. El nuevo borde esta part_in_clu posiciones a la derecha
       ! del anterior
       bordes_clusters(iclus)=bordes_clusters(iclus-1)+part_in_clu
    enddo ! Outmost loop: Cluster nuevo

    !         write(11,997)ijij,ksi,aux,aux ! configuracion, nro. de clusters, basura, basura

    ! Calculo cantidad total de fragmentos (se usa como interruptor)
    nclus0=sum(mass_dist)
    return
  end subroutine mst
  !______________________________________________________________________________

  function epot(d2,isos1,isos2)
    Real :: d2,epot
    logical :: isos1,isos2
    Real    :: rij,eint

    rij=sqrt(d2)
    eint=0.0
    if(isos1.neqv.isos2) then      ! vnp
       eint=(vr*exp((-mur)*rij))/rij-va*(exp((-mua)*rij))/rij-ep0_np
    else
       eint=(v0*exp((-mu0)*rij))/rij-ep0_nn
    endif
    if(isos1.eqv..true. .and. isos2.eqv..true.) eint=eint+ee*(exp(-invlambda*rij))/rij
    epot=eint

    return
  end function epot


  !______________________________________________________________________________
end module mod_mst
