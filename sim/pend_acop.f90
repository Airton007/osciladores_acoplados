program osc_acopladas_n

    ! Este programa simula o movimento de um sistema de N pêndulos acoplados, resolvendo numericamente 
    ! as equações diferenciais que descrevem o sistema usando o método de Runge-Kutta de quarta ordem.
    ! O programa permite que o usuário escolha entre três conjuntos de condições iniciais diferentes e 
    ! gera um arquivo de saída com os resultados da simulação.
    
      use iso_fortran_env, only: i4 => int32, i8 => int64, sp => real32, dp => real64
      use runge_kutta ! Importa o módulo contendo a implementação do método de Runge-Kutta
      implicit none
    
      ! Declaração de variáveis
      real(kind=dp) :: t, t0, h, omega0         ! Parametros de discretização e frequência natural
      real(kind=dp), allocatable :: z(:)        ! Vetor de estado do sistema
      integer(kind=i8) :: N, i, j, arquivo      ! Variaveis auxiliares
      integer(kind=i4) :: d, N_pendulos         ! Dimensão do sistema e número de péndulos
      character(len=50) :: nome_arquivo
    
      ! Parâmetros do problema
      real(kind=dp), parameter :: m = 1.0_dp, k = 1.0_dp  ! Massa e constante elástica
      real(kind=dp), parameter :: g = 9.81_dp, l = 1.0_dp  ! Gravidade e comprimento da corda
    
      ! Frequência natural do sistema
      omega0 = sqrt(g / l)
    
      ! Entrada do número de pêndulos
      print*, 'Digite o número de pêndulos acoplados (N):'
      read*, N_pendulos
    
      if (N_pendulos < 2) then
         print*, 'O número de pêndulos deve ser pelo menos 2.'
         stop
      end if
    
      ! Definições para a integração
      d = 2 * N_pendulos           ! Dimensão do sistema (posição e velocidade para cada pêndulo)
      N = 10000                    ! Número de passos
      h = 0.001_dp                 ! Passo de tempo
      t0 = 0.0_dp                  ! Tempo inicial
      t = t0
    
      ! Alocação do vetor z
      allocate(z(d))
    
      ! Escolha das condições iniciais
      print*, 'Escolha as condições iniciais:'
      print*, '1. Todas as posições iguais'
      print*, '2. Posições alternadas (opostas)'
      print*, '3. Posições diferentes'
      read*, i
    
      if (i == 1) then
         ! Condições iniciais: todas as posições iguais
         z(1:N_pendulos) = 1.0_dp       ! Posições iniciais
         z(N_pendulos+1:d) = 0.0_dp     ! Velocidades iniciais
         nome_arquivo = 'pend_acop_pos_iguais.dat'

      elseif (i == 3) then
         ! Condições iniciais: posições diferentes
         z(1:N_pendulos) = [(1.0_dp * j, j=1, N_pendulos)]
         z(N_pendulos+1:d) = 0.0_dp
         nome_arquivo = 'pend_acop_pos_diferentes.dat'
      
      elseif (i == 3) then
         ! Condições iniciais: posições alternadas (opostas)
         z(1:N_pendulos) = [((-1.0_dp)**j, j=1, N_pendulos)]
         z(N_pendulos+1:d) = 0.0_dp
         nome_arquivo = 'pend_acop_pos_opostas.dat'
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
          write(arquivo, *) t, z
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
    
        integer(kind=i4) :: j
    
        ! Derivadas das posições
        ! A derivada da posição (dx/dt) é igual à velocidade. 
        ! Aqui, a velocidade é armazenada na parte superior do vetor z, do índice N_pendulos+1 até d.
        res(1:N_pendulos) = z(N_pendulos+1:d)  ! dx/dt = v

        ! Derivadas das velocidades
        ! Calcula a derivada da velocidade (dv/dt) para cada pêndulo, considerando as forças atuantes.
        do j = 1, N_pendulos
           if (j == 1) then
              ! Primeiro pêndulo: 
              ! A força restauradora inclui o termo da mola com o próximo pêndulo (z(j+1) - z(j)).
              res(N_pendulos+j) = (-omega0**2 * z(j) + k * (z(j+1) - z(j))) / m
           elseif (j == N_pendulos) then
              ! Último pêndulo: 
              ! A força restauradora inclui o termo da mola com o pêndulo anterior (z(j-1) - z(j)).
              res(N_pendulos+j) = (-omega0**2 * z(j) + k * (z(j-1) - z(j))) / m
           else
              ! Pêndulos intermediários:
              ! A força restauradora inclui os termos das molas com os pêndulos adjacentes
              ! (z(j-1) + z(j+1) - 2 * z(j)) representa a interação com os vizinhos.
              res(N_pendulos+j) = (-omega0**2 * z(j) + k * (z(j-1) + z(j+1) - 2.0_dp * z(j))) / m
           end if
        end do

      end function
    
    end program
    