
program osc_acopladas

! Este programa simula o movimento de um sistema de dois pêndulos acoplados, resolvendo numericamente 
! as equações diferenciais que descrevem o sistema usando o método de Runge-Kutta de quarta ordem. 
! O programa permite que o usuário escolha entre três conjuntos de condições iniciais diferentes e 
! gera um arquivo de saída com os resultados da simulação.

    use iso_fortran_env, only: i4 => int32, i8 => int64, sp => real32, dp => real64
    use runge_kutta ! Importa o módulo contendo a implementação do método de Runge-Kutta
    implicit none
  
    ! Declaração de variáveis
    real(kind=dp) :: t, t0, h, omega0       ! Parametros de discretização e frequência natural
    real(kind=dp), allocatable :: z(:)      ! Vetor de estado do sistema
    integer(kind=i8) :: N, i, arquivo       ! Variaveis auxiliares
    integer(kind=i4), parameter :: d = 4    ! Dimensão do sistema
    character(len=50) :: nome_arquivo
    
    ! Parâmetros do problema
    real(kind=dp), parameter :: m = 1.0_dp, k = 1.0_dp  ! Massa e constante elástica
    real(kind=dp), parameter :: g = 9.81_dp, l = 1.0_dp  ! Gravidade e comprimento da corda
  
    ! Frequência natural do sistema
    omega0 = sqrt(g / l)
  
    ! Definições para a integração
    N = 1000000                 ! Número de passos
    h = 0.0001_dp               ! Passo de tempo
    t0 = 0.0_dp                 ! Tempo inicial
    t = t0
  
    ! Alocação do vetor z
    allocate(z(d))
  
    ! Escolha das condições iniciais
    print*, 'Escolha as condições iniciais:'
    print*, '1. x1 = x2'
    print*, '2. x1 = -x2'
    print*, '3. x1 diferente de x2'
    read*, i
  
    if (i == 1) then
       ! Condições iniciais para x1 = x2
       z(1) = 1.0_dp     ! x1 posição inicial
       z(2) = 1.0_dp     ! x2 posição inicial
       z(3) = 0.0_dp     ! v1 velocidade inicial
       z(4) = 0.0_dp     ! v2 velocidade inicial
       nome_arquivo = 'pend_acop_x1_eq_x2.dat'
    elseif (i == 2) then
       ! Condições iniciais para x1 = -x2
       z(1) = 1.0_dp     
       z(2) = -1.0_dp   
       z(3) = 0.0_dp     
       z(4) = 0.0_dp     
       nome_arquivo = 'pend_acop_x1_eq_neg_x2.dat'
    elseif (i == 3) then
       ! Condições iniciais para x1 diferente de x2
       z(1) = 1.0_dp     
       z(2) = 0.0_dp     
       z(3) = 0.0_dp    
       z(4) = 0.0_dp     
       nome_arquivo = 'pend_acop_x1_diff_x2.dat'
    else
       print*, 'Opção inválida.'
       stop
    end if
  
    ! Abertura de arquivo para salvar resultados
    open(newunit=arquivo, file=nome_arquivo, status="replace")
  
    ! Loop principal para integração numérica usando Runge-Kutta
    do i = 0, N-1
        t = t0 + i * h                           
  
        ! Atualiza as variáveis com Runge-Kutta
        z = rk4_multi(t, z, h, d, derivada_z)        
        write(arquivo, *) t, z(1), z(2), z(3), z(4) 
    end do
    close(arquivo)
    deallocate(z)
  
  contains
  
    ! Função para calcular as derivadas do sistema
    function derivada_z(t, z, d) result(res)
      ! Parâmetros locais da função
      real(kind=dp), intent(in) :: t
      real(kind=dp), intent(in) :: z(:) 
      integer(kind=i4), intent(in) :: d
      real(kind=dp) :: res(d)
  
      ! Derivadas das posições
      res(1) = z(3)  ! dx1/dt = v1
      res(2) = z(4)  ! dx2/dt = v2
  
      ! Derivadas das velocidades
      res(3) = (-omega0**2 * z(1) + k * (z(2) - z(1))) / m
      res(4) = (-omega0**2 * z(2) + k * (z(1) - z(2))) / m
    end function
  
  end program