% Carregar a toolbox
addpath('G:\\FDAfuns');  % Substitua pelo caminho real da pasta 'fdasrvf_MATLAB-master'
% Consulte a biblioteca para dúvidas

% Dados
Y = ; % Seus dados
fs = 600;
dt = 1/fs;
t =  dt:dt:size(Y, 2)/fs;

% Configurar a base B-spline
%Número de bases (200): Este é o número de funções de base usadas para criar a representação suavizada da curva. Ao ajustar os dados a um modelo de função de base, estamos basicamente tentando representar nossos dados como uma combinação de funções de base (como polinômios, splines, etc.). O número de funções de base define o número máximo de características distintas que o modelo pode capturar nos dados. Escolher o número certo de bases é essencial: muito poucas podem não capturar adequadamente a complexidade dos dados, enquanto muitas podem resultar em overfitting.
%Ordem das bases (4): Este é o grau das funções de base usadas. Por exemplo, se você está usando funções de base polinomiais, uma ordem de 4 significaria que a função de base mais complexa usada seria um polinômio de quarto grau.
%Mínimo de número harmônicos (2): Este é o número mínimo de harmônicos usados na análise. Na análise harmônica, que é um tipo de análise de séries temporais, os harmônicos são funções sinusoidais usadas para decompor a série temporal em componentes de frequência diferentes.
%Taxa de amostragem (600 Hz): Este é o número de amostras coletadas por segundo. Uma taxa de amostragem mais alta pode capturar mais detalhes do sinal, mas também pode levar a mais ruído e requer mais recursos computacionais.
%Lambda (1e-6): Este é um parâmetro de suavização. A suavização é usada para reduzir o ruído nos dados. Um lambda maior levará a uma suavização maior (isto é, uma curva mais "suave"), enquanto um lambda menor resultará em menos suavização. A escolha do lambda pode ser um equilíbrio entre reduzir o ruído e preservar os detalhes dos dados.

nBasis = 200; % Ajuste este valor com base na complexidade dos dados
order = 4;
t_range = [min(t), max(t)]; %necessário para gerar o tamanho do sinal suavizado
basis = create_bspline_basis(t_range, nBasis, order); % necessário para gerar o objeto fd

% Criar objeto fd vazio
coef = zeros(nBasis, size(Y, 1));
fdobj = fd(coef, basis);            %base para gerar o sinal suavizado

% Suavização e conversão para objeto funcional
lambda = 1e-6; % Parâmetro de regularização para suavização
Lfdobj = int2Lfd(2); % Parâmetro de transformação do dado
fdParobj = fdPar(fdobj, Lfdobj, lambda);
smooth_data = smooth_basis(t, Y', fdParobj); % dados suavizados para extrair os sinais correspondentes aos eventos

% FPCA
nharm = ; % número de harmônicos para extrair com base na complexidade de eventos da tarefa
pca_results = pca_fd(smooth_data, nharm);   % análise de componentes principais ajustada para dados funcionais


% Plotar a média, variância e autofunções
figure;
subplot(2, 2, 1);
plot(t, eval_fd(t, pca_results.meanfd));
title('Média');
xlabel('Tempo');
ylabel('Valor');

subplot(2, 2, 2);
plot(pca_results.values);
title('Valores próprios');
xlabel('Número da componente principal');
ylabel('Valor próprio');
xlim ([1 4]);

for i = 1:nharm
    subplot(2, 2, i + 2);
    plot(t, eval_fd(t, pca_results.harmfd(i)));
    title(['Autofunção ' num2str(i)]);
    xlim ([10 590]);
    ylim ([-5e-1 5e-1]);
    xlabel('Tempo');
    ylabel('Valor');
end
