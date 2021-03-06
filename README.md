# Trabalho 1
Um algoritmo exato branch-and-bound e heurísticas para o problema da mochila múltipla

## Implementação
As seguintes tarefas de programação devem ser realizadas, usando as rotinas da biblioteca GLPK:
- [x] I1: implemente um algoritmo de branch-and-bound para o MKP usando a formulação (F1);
- [x] I2: Proponha e implemente pelo menos duas heurísticas ingênuas para gerar boas soluções para o
MKP. Dica: discutimos, na aula prática, duas heurísticas simples, uma gulosa e outra aleatória, para gerar soluções para o MKP. A ideia é propor alguma melhoria em alguma dessas heurísticas, ou propor uma outra heurística simples (sem o uso de metaheurísticas)

## Testes
- [x] T1: (30% da nota) Verifique a qualidade das soluções geradas pela implementação I1, limitado a um
tempo limite de 5 minutos. Para medir a qualidade da solução, calcule o gap de dualidade, que é
dado por 100%(UB − LB)/UB e dá uma garantia para a qualidade da solução. Quando o gap é
zero, temos que a solução encontrada é ótima. Reporte o gap médio, o total de instâncias de testes
resolvidas na otimalidade e o tempo médio gasto nesses casos;
- [x] T2: (10% da nota) Verifique a qualidade dos limitantes superiores dados pela relaxação linear, usando a
implementação I1. Para medir a qualidade dos limitantes UB obtidos pela relaxação linear, calcule
o gap de otimalidade, que é dado por 100%(UB − z∗)/z∗ e dá uma medida para a qualidade do
limitante considerando-se o valor da melhor solução conhecida como o valor de z∗. O valor de
z∗ a ser usado nesta análise deve ser o valor da melhor solução gerada pela implementação I1 no
experimento T1. Quanto mais próximo de zero o gap está, melhor é a qualidade do limitante gerado
pela relaxação linear. Reporte o gap médio, o total de instâncias de testes em que o gap foi zero e
o tempo médio gasto pela relaxação linear em todas as instâncias testes;
- [x] T3: (60% da nota) Verifique a qualidade das soluções geradas pelas heurísticas propostas e implementadas em I2. Calcule a qualidade das soluções, usando a fórmula do gap de dualidade e considerando
o valor da relaxação linear como UB da instância. Reporte o gap médio de cada heurística, o total
de instâncias testes em que cada heurística obteve melhor resultado dentre as heurísticas propostas
e o tempo médio gasto pela heurística.

# Trabalho 2

## Implementação
- [x] I1: implemente um algoritmo de branch-and-bound para o MKP usando a formulação (F1);
- [x] I2: Proponha e implemente uma heurística baseada em modelo matemático da classe das heurística de
melhoria ou da classe de heurística baseadas na relaxação para gerar boas soluções para o MKP.
Dica: discutimos, durante as aulas, algumas heurísticas dessas classes para o MKP. A
ideia é implementar alguma dessas heurísticas, ou propor uma outra heurística usando
um dos métodos dessa classe visto em aula.

## Testes

- [x] T1: (80% da nota) Verifique a qualidade das soluções geradas pelas heurísticas propostas e implementadas em I2. Calcule a qualidade das soluções, usando a fórmula do gap de dualidade, que é computado
pela fórmula 100(UB − LB)/UB e dá uma garantia para a qualidade da solução. Quando o gap
é zero, temos que a solução encontrada é ótima. Nessa fórmula, LB é o valor da solução gerada
pela heurística e UB é o valor da relaxação linear (ou o valor do melhor limitante dual obtido pelo
algoritmo exato na implementação I1). Reporte o gap médio de cada heurística, o total de instâncias
testes em que cada heurística obteve melhor resultado dentre as heurísticas propostas e o tempo
médio gasto pela heurística.
- [x] T2: (20% da nota) Compare os resultados de desempenho e de qualidade das soluções geradas pela
implementação I2 em relação à qualidade e desempenho das heurísticas ingênuas implementadas
no primeiro trabalho prático. Caso o grupo não tenha implementado as heurísticas ingênuas, fica
a critério do grupo implementá-las nesta etada ou usar o código das heurísticas ingênuas a ser
disponibilizado pelo professor no Moodle.
