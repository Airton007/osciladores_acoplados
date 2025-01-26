module runge_kutta
    ! Este módulo implementa o método de Runge-Kutta de quarta ordem (RK4), uma técnica popular para resolver sistemas de equações diferenciais.
    ! O RK4 é bastante utilizado por sua precisão e eficiência, sendo ideal para simulações de sistemas dinâmicos.
    ! O módulo contém uma interface para a função de derivadas, que deve ser fornecida pelo usuário, e uma função principal (rk4_multi) 
    ! que realiza a integração numérica do sistema, passo a passo.

    use iso_fortran_env, only : i4 => int32, i8 => int64, sp => real32, dp => real64
    implicit none

    ! Aqui, definimos uma interface para a função de derivadas, que calcula as taxas de variação das variáveis do sistema.
    ! Essa função será fornecida pelo usuário e deve seguir o formato esperado para funcionar corretamente.
    interface
        function derivada_multi(t, x, d) result(res)
            ! Utilizamos os tipos padrão do Fortran para garantir consistência no código.
            use iso_fortran_env, only : i4 => int32, i8 => int64, sp => real32, dp => real64
            implicit none 
            integer(i4), intent(in) :: d   ! Número de variáveis ou dimensões do sistema.
            real(dp), intent(in) :: t     ! O tempo no momento atual.
            real(dp), intent(in) :: x(:)  ! Vetor contendo as variáveis de estado do sistema.
            real(dp) :: res(d)            ! Resultado das derivadas, um vetor com "d" elementos.
        end function
    end interface
    
contains

    ! A função rk4_multi é a principal do módulo e realiza a integração do sistema de equações diferenciais.
    ! Ela aplica o método de Runge-Kutta de quarta ordem para calcular a evolução do sistema a cada passo de tempo.
    function rk4_multi(ti, xi, h, d, g) result(res)
        procedure(derivada_multi) :: g  ! A função "g" calcula as derivadas e é fornecida pelo usuário.
        integer(i4), intent(in) :: d    ! O número de dimensões do sistema (quantidade de variáveis a serem integradas).
        real(kind=dp) :: ti             ! O tempo inicial, ou seja, o ponto de partida da integração.
        real(kind=dp) :: h              ! O passo de integração, que determina o tamanho de cada "passo" no tempo.
        real(kind=dp) :: xi(d)          ! O estado inicial do sistema (valores das variáveis no tempo ti).
        real(kind=dp) :: k1(d), k2(d), k3(d), k4(d)  ! Os coeficientes intermediários usados pelo método RK4.
        real(kind=dp) :: res(d)         ! O resultado final, ou seja, o novo estado do sistema após o passo de integração.

        ! Cálculo dos coeficientes intermediários k1, k2, k3 e k4, que são usados para estimar a evolução do sistema:
        k1 = h * g(ti, xi, d)                              ! Primeiro coeficiente, avaliado no estado inicial.
        k2 = h * g(ti + h/2.0_dp, xi + k1/2.0_dp, d)       ! Segundo coeficiente, calculado com o estado ajustado por k1.
        k3 = h * g(ti + h/2.0_dp, xi + k2/2.0_dp, d)       ! Terceiro coeficiente, ajustado por k2.
        k4 = h * g(ti + h, xi + k3, d)                     ! Quarto coeficiente, calculado com o estado ajustado por k3.

        ! Combinamos os coeficientes para obter a estimativa do próximo estado do sistema.
        res = xi + (k1 + 2.0_dp*k2 + 2.0_dp*k3 + k4)/6.0_dp
    end function

end module runge_kutta
