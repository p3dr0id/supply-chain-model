#Definição das Variáveis
using DifferentialEquations
using LinearAlgebra
num_elos = 2 #Número de elos da cadeia de suprimentos
tempo_adap=1 #Tempo de adaptação
tau=0.01292929 #Tempo de antecipação de estoque
epsilon=8.990909 #Tempo de ajuste de estoque
beta=1 #Tempo de ajuste de produção
fornecedor = (1:1:num_elos) #cria vetor de 1 ao n° de osciladores com passo 1

num_sim=1.0 #Define número médio de simulações para cálculo da média das medidas

tf = 4.00 #Define tempo final de simulação
dt=1 #Define intervalo entre pontos de tempo de simulação
tspan=(0.0,tf)#Define vetor de tempo de simulação
#println(typeof(tspan))
#println(tspan)

n0=ones(1,num_elos) #Define condição inicial para os níveis de estoque
q0=zeros(1,num_elos) #Define condição inicial para as taxas de entrega
v0=zeros(1,num_elos) #Define condição inicial para as derivadas das taxas de entrega

R = [zeros(num_elos-1,1);[-2.0];zeros(num_elos,1); zeros(num_elos-1,1);[1.0]] #vetor coluna de zeros com dimensao No sendo o último elemento da coluna 1 (associado ao consumidor)
#R=vcat(zeros(num_elos-1,1),[-2.0],zeros(num_elos,1), zeros(num_elos-1,1),[1.0])
println("R:")
println(typeof(R))
println(R)

variabilidade_sim = Float64[]
gamma = ((beta + epsilon)/tempo_adap) / 2
omega2 = 1/(tempo_adap * tau)

#definição da função de nível de estoque e taxa de entrega para os elos da cadeia
function f(theta,p, t)
    N = num_elos
    #N = floor(size(theta) / 3)
    n = theta[1:N]
    q = theta[N+1:2*N]
    v = theta[2*N+1:3*N]
    dndt = zeros(N)
    f = zeros(N)
    a = zeros(N)
    for i in range(1,N-1,step=1)
        dndt[i] = q[i] - q[i + 1]
        f[i]=(1/tempo_adap)*((q[i+1]/tau)+beta*v[i+1])
        a[i]= f[i]-2*gamma*v[i]-omega2*q[i]
    end
    a[N]= -2*gamma*v[N]-omega2*q[N]
    #return [dndt, v, a]
    return vcat(dndt,v,a)
end

#definição da função que acrescenta aleatoriedade ao sistema
function g(theta,p,t)
    algumacoisainteressante=0.1
    cons_estoc=R*algumacoisainteressante
    return cons_estoc
end

yinit=hcat(n0, q0, v0)
#yinit=[ 1. 1. 0. 0. 0. 0.]
println("Y init:")
println(typeof(yinit))
println(yinit)

#for sim in 1:num_sim
#integração do sistema de equações diferenciais
prob = SDEProblem(f,g,yinit,tspan)
sol = solve(prob,EM(),dt=dt)
println(sol)
# using Plots
# p=plot(sol)
# savefig(p,"cadeia.png")
