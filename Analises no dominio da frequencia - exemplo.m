%% CARREGA OS DADOS DE LFP
clear all
close all
clc

% Exemplo 
% Carregue o banco de dados
   addpath('C:\Registros') 
   %exemplo 'load REGISTRO_LAGARTOS'
   load X % exemplo para núcleo X
   load Y % exemplo para núcleo Y
% Carregue as funções e scripts necessários para limpar os dados 
% (que não seja nativo do MATLAB)

    %Importe de https://github.com/tortlab/SignalAnalysis2020.2
    addpath('C:\SignalAnalysis2020.2');
    
    %Para eegfilt.m - EEGLAB: http://sccn.ucsd.edu/eeglab/
    addpath('C:\eeglab')
    
    
% Os exemplos foram criados com base em um núcleo, caso necessário substituir X por Y    
%% Visualize e remova os registros ruidosos
a = 'Canal'; % código do canal a ser inspecionado
srate = 600; % taxa de amostragem de coleta
idx = 1/srate; % passo de incremento no tempo
data = X(a,:); % dado selecionado a partir do canal
timevector = idx:idx:length(data)/srate; %vetor de tempo gerado do valor inicial até o valor final em passos de 'idx'
plot(timevector,data); %plot do sinal para investigar os ruídos

%% Remova os ruídos se forem necessários ou os declare
% Canais registrados e áreas 
% Declare o canal a sua respectiva área de posicionamento
% 1, 2, 3, 4 ,5 ...

X_SR = X(1,[1:t1*srate,t2*600:srate]); %X_SR será o vetor X sem ruído (X_SR)

%% Filtros de 60 Hz e menores que 1 Hz
% Bandstop filter de pulso infinito
% Filtro de 60 Hz
d = designfilt('bandstopiir','FilterOrder',10, ...
    'HalfPowerFrequency1',59,'HalfPowerFrequency2',61, ...
    'DesignMethod','butter','SampleRate',600); 
% consulte a descrição de designfilt para ajustes na função

X_F = filtfilt(d,X_SR(1,:)); % Sinal filtrado a partir do vetor X_SR

% Filtro para bandas <1 Hz   
clear d % limpar o filtro 'd'
d = designfilt('highpassiir', 'FilterOrder',5, ...
    'PassbandFrequency',0.9, ...
    'PassbandRipple',0.5,'SampleRate',600);
% consulte a descrição de designfilt para ajustes na função

X_F = filtfilt(d,X_F(1,:)); % Sinal filtrado a partir do vetor X_F
% Y_F = filtfilt(d,Y_F(1,:)); % Sinal filtrado caso queira parear X e Y -
% Y_F é hipotético
%% Declarar bandas de frequência
DATA = X_F(1,:);
srate = 600; % frequencia de amostragem em Hz
dt = 1/srate; % Passo em segundos
time_vector = dt:dt:length(DATA)/srate; % Vetor de tempo em segundos
LFP = DATA;

% Decomposição da frequência no tempo 
frequency_vector = 1:0.1:100;
frequency_width = 1;
clear TFDcont
TFDcont = zeros(length(frequency_vector),length(LFP));

% Loop para criar vetor de frequêncais ao longo do tempo
for j=1:length(frequency_vector)
    filtrado = eegfilt(LFP,srate,frequency_vector(j),...
        frequency_vector(j)+frequency_width);
    
    TFDcont(j,:)= abs(hilbert(filtrado ));
end

%FILTRA O SINAL a partir da função eegfilt
delta = eegfilt(LFP,srate,1,4);
theta = eegfilt(LFP,srate,6,12);
beta  = eegfilt(LFP,srate,15,25);
gamma = eegfilt(LFP,srate,25,55);

%OBTEM O ENVELOPE DE AMPLITUDE a partir da função nativa
% abs (valores absolutos - apenas positivos)
% a transformada de hilbert declara valores discretos no sinal contínuo e
% extrai valores de amplitude e frequência
deltaAmp = abs(hilbert(delta));
thetaAmp = abs(hilbert(theta));
betaAmp = abs(hilbert(beta));
gammaAmp = abs(hilbert(gamma));

%Agrupamos e nomeamos as bandas
Bands={'Delta','Teta','Beta','Gama'};

% Gerar figura dos dados
close all
fig1 = figure(1);clf
set(gcf,'color','white')
subplot(2,2,1)
plot(time_vector,LFP,'k','linewidth',1)
set(gca,'fontsize',12)
ylim([-1 1])
x11 = gca;
ylabel('Voltage (mV)')
xlabel('Tempo (s)')
box off

subplot(2,2,3)
imagesc(time_vector,frequency_vector+frequency_width/2,(TFDcont))
axis xy
xlabel('Tempo (s)')
ylabel('Frequência (Hz)')
xlim([1 600])
ylim([1 50])
caxis([0 0.05])
set(gca,'fontsize',12)
box off
title('Decomposição de frequências no tempo','fontsize',16)
x12 = gca;

subplot(2,2,[2 4])
var=0.2;
plot(time_vector,delta+6*var,'b','linewidth',1)
hold on
plot(time_vector,deltaAmp+6*var,'b','linewidth',2)

plot(time_vector,theta+3*var,'c-','linewidth',1)
plot(time_vector,thetaAmp+3*var,'c-','linewidth',2)

plot(time_vector,beta-3*var,'g-','linewidth',1)
plot(time_vector,betaAmp-3*var,'g-','linewidth',2)

plot(time_vector,gamma-6*var,'m-','linewidth',1)
plot(time_vector,gammaAmp-6*var,'m-','linewidth',2)

plot(time_vector,LFP-9*var,'k-','linewidth',1)
plot(time_vector,LFP-9*var,'k-','linewidth',1)

title('Decomposição dos sinais','fontsize',16)
set(gca,'fontsize',12)
xlim([295 298])
ylim([-2.4 2.4])
xlabel('Tempo(s)')
set(gca,'YTick',[-6*var:0.4:6*var],'YTickLabel',flip(Bands))
linkaxes([x11 x12],'x');
box off


%% Frequência relativa, Espectro de densidade e Espectrograma
clc
clear  Psd*  T* ACG* P* Spect* freq t

    clear data* PSDT
    data   = X_F(1,:);       % dados coletado para análise

%% Frequência relativa
% Definição das bandas de frequência mais presentes
faixa_delta = [1 4];
faixa_theta = [6 12];
faixa_alpha = [7 12];
faixa_beta  = [15 25];
faixa_gama1  = [25 55];

% Encontra os índices correspondentes às faixas de frequência
idx_delta = find(F >= faixa_delta(1) & F <= faixa_delta(2));
idx_theta = find(F >= faixa_theta(1) & F <= faixa_theta(2));
idx_alpha = find(F >= faixa_alpha(1) & F <= faixa_alpha(2));
idx_beta = find(F >= faixa_beta(1) & F <= faixa_beta(2));
idx_gama1 = find(F >= faixa_gama1(1) & F <= faixa_gama1(2));

% Extrai os valores de frequência relativa para cada banda de frequência
clear frequencia_*
for i = 1:5
frequencia_delta(i,:) = data(i,idx_delta);
frequencia_theta(i,:) = data(i,idx_theta);
frequencia_alpha(i,:) = data(i,idx_alpha);
frequencia_beta(i,:) =  data(i,idx_beta);
frequencia_gama1(i,:) = data(i,idx_gama1);

end
% Calcula a média dos canais para cada banda de frequência
media_delta = mean(mean(frequencia_delta, 2));
media_theta = mean(mean(frequencia_theta, 2));
media_alpha = mean(mean(frequencia_alpha, 2));
media_beta = mean(mean(frequencia_beta, 2));
media_gama1 = mean(mean(frequencia_gama1, 1));

% Normaliza e define a proporção das bandas de frequência
delta = media_delta/sum([media_delta,media_theta,media_beta,media_gama1]);
teta = media_theta/sum([media_delta,media_theta,media_beta,media_gama1]);
beta = media_beta/sum([media_delta,media_theta,media_beta,media_gama1]);
gama = media_gama1/sum([media_delta,media_theta,media_beta,media_gama1]);

% Plot do gráfico das médias para cada banda de frequência
bandas = {'Delta', 'Teta', 'Beta', 'Gama'};
medias = [delta; teta; beta; gama];
figure()
subplot(2,1,1)
bar(medias);
set(gca, 'XTickLabel', bandas);
xlabel('Banda de Frequência');
ylabel('`Proporção relativa das bandas');
title('Média do CDM para Cada Banda de Frequência');

%% Espectro de densidade e espectrograma
clc
clear  P* T* idx F data*

data       = (X_F(1,:)); %tempo total
data1      = (X_F(1,1:length(data)/2));             % Etapa de ambientação (AMB)
data2      = (X_F(1,length(data)/2:length(data)));  % Etapa de exposição   (EXP)
srate      = 600;        % taxa de amostragem / janela de tempo
idx        = 1/srate;    % incremento de tempo
window     = 5*srate;    % tamanho da janela de análise, é interessante ajustar a janela a partir do 'srate'
overlap    = window*0.2; % sobreposição feita a partir da próxima janela (window), para não perder a continuidade
nfft       = 2^13;       % número de pontos que serão coletados a partir da transformada de fourier

%AMB
% Densidade de espectro de potência da ambientação
[PSD1,F]   = pwelch(data1,window,overlap,nfft,srate);        %função de welch que estrai valores de densidade de potência das frequências
% Espectrograma da ambientação
[S1 F T1 P2]  = spectrogram(data,window,overlap,nfft,srate); %função que extrai valores a partir da transformada curta de fourier

%EXP
% Densidade de espectro de potência da exposição
[PSD2,F]   = pwelch(data2,window,overlap,nfft,srate); 
% Espectrograma da exposição
[S2 F T2 P2]  = spectrogram(data,window,overlap,nfft,srate);

%% Análise dos valores de potência entre as condições AMB vs. EXP
%  TESTE T NOS DADOS NORMALIZADOS PELA MEDIA - TEMPO TOTAL
clc
% Bands=[2,6; 7,12; 15,25; 36,40];%definir bandas de interesse
Bands=[2,4; 7,12; 15,25; 29,33; 36,40];%definir bandas de interesse
clear h* p s IC NormAMB NormEXP idx TGH
for b = 1:5 % número de bandas
    clear idx Norm
    idx          = find(F>Bands(b,1) & F<Bands(b,2));
    Norm         = mean([mean(PSD1(:,idx),2), mean(PSD2(:,idx),2)],2); %fator de nomalização - média das médias de AMB e EXP
    NormAMB(:,b) = mean((PSD1(:,idx)),2)./Norm; % AMB normalizado
    NormEXP(:,b) = mean((PSD2(:,idx)),2)./Norm; % EXP normalizado   
    TGH{b} = ([NormAMB(:,b);NormEXP(:,b)]);  % Une os dados para investigar normalidade e homocedacidade
    [hv{b},pv{b},adstat{b},cv{b}] = adtest(TGH{b}); % Investiga a normalidade
    group = ([1;1;1;1;1;2;2;2;2;2]); % vetor de fatores AMB = 1 e EXP = 2
    vartestn(TGH{b},group,'TestType','Bartlett'); % Investiga a homocedacidade
    [h(b) p(b) ci stats]  = ttest(NormAMB(:,b), NormEXP(:,b),'Alpha',0.05,'Tail','left');
    IC(:,:,b) = ci; % Intervalo de confiança
    s{b} = stats;   % Valores da estatística
    effect(:,b) = computeCohen_d(NormAMB(:,b), NormEXP(:,b)); % Tamanho do efeito 'Cohen'
    
    % Investiga se a banda 01 possui potência maior do que as bandas
    % subsequentes
    [h2(b) p2(b) ci2 stats2]  = ttest(NormAMB(:,1), NormAMB(:,b),'Alpha',0.05,'Tail','both');
    IC2(:,:,b) = ci2;
    s2{b} = stats2;
    effect2(:,b) = computeCohen_d(NormAMB(:,1), NormAMB(:,b), 'paired');
    
end
p([1,2,3,4,5]) % Apresenta os valores p para cada uma das bandas

%% Preparar o plot dos dados com base na densidade (log), no espectro, na decomposiçao no tempo e nos valores de potência de AMB e EXP
% Espectrograma (Potencia no tempo)
SpectAMB(:,:,a) = P1(:,1:T1); %AMB
SpectEXP(:,:,a) = P2(:,1:T2); %EXP

% Decomposição de frequência no tempo
frequency_vector = 1:0.1:59; % vetor de frequências
frequency_width = 1;         % valor de incremento para o sinal filtrado
clear TFDcont
TFDcont = zeros(length(frequency_vector),length(data(1,:)));

% Loop para criar um vetor de frequências ao longo do tempo para plot
for j=1:length(frequency_vector)
    filtrado = eegfilt(mean(data(:,:),1),srate,frequency_vector(j),...
        frequency_vector(j)+frequency_width);
    TFDcont(j,:)= abs(hilbert(filtrado)); % vetor para o plot de frequências
end

%% Figura

FIG1=figure(1);clf
subplot(4,4,[1 4])
LFP_1 = data1; %sinal filtrado AMB
LFP_2 = data2; %sinal filtrado EXP
dt = 1/srate;
time_vector = dt:dt:length(LFP_1)/srate;
plot(time_vector,LFP_1,'k','linewidth',1)   %sinal no tempo para AMB
hold on                                     %sobreposição
plot(time_vector,LFP_2-1,'r','linewidth',1) %sinal  no tempo para EXP
ylim([-1.8 .5])
box off
xlabel 'Tempo (s)'
x1=195;
xlim([x1 x1+3])
set(gca,'ytick',[-1 0],'yticklabels',{'Exposição','Ambientação'})
set(gca, 'Fontsize', 10)

%Filtro seletivo das bandas
clear delta* alfa* beta* gama*
delta_AMB = eegfilt(LFP_1,srate,1,4);
delta_EXP = eegfilt(LFP_2,srate,1,4);

theta_AMB = eegfilt(LFP_1,srate,6,12);
theta_EXP = eegfilt(LFP_2,srate,6,12);

beta_AMB = eegfilt(LFP_1,srate,12,25);
beta_EXP = eegfilt(LFP_2,srate,12,25);

gama_AMB = eegfilt(LFP_1,srate,30,50);
gama_EXP = eegfilt(LFP_2,srate,30,50);

subplot(4,4,[1 4]) %continuando a sobreposição no mesmo conjunto do sinal bruto
plot(time_vector,delta_AMB-0.2,'k','linewidth',1)
plot(time_vector,delta_EXP-1.2,'r','linewidth',1)
plot(time_vector,theta_AMB-0.35,'k','linewidth',1)
plot(time_vector,theta_EXP-1.35,'r','linewidth',1)
plot(time_vector,beta_AMB-0.5,'k','linewidth',1)
plot(time_vector,beta_EXP-1.5,'r','linewidth',1)
plot(time_vector,gama_AMB-0.5,'k','linewidth',1)
plot(time_vector,gama_EXP-1.5,'r','linewidth',1)

% Espectrograms
subplot(4,4,[5 6 7])
imagesc(TAMB{1},F,mean(SpectAMB,3))              % Espectrograma da AMB
hold on
imagesc(TEXP{1}+TAMB{1}(end),F,mean(SpectEXP,3)) % Espectrograma da EXP
axis xy
ylim([20 40])                                    % Eixo das frequências
xlim([0 600])                                    % Eixo do tempo
caxis([0 0.005])                                 % Eixo do valor de potência (ajustar)
title 'CDM - AMB.'
xlabel 'Tempo (s)'
ylabel 'Frequência (Hz)'
colorbar
set(gca, 'Fontsize', 10)

subplot(4,4,[9 10 11])
imagesc(TAMB{1},F,mean(SpectAMB,3))               % Espectrograma de AMB
hold on
imagesc(TEXP{1}+TAMB{1}(end),F,mean(SpectEXP,3))  % Espectrograma de EXP
axis xy
ylim([1 20])                                      % Eixo Y das frequências
xlim([0 600])                                     % Eixo X do tempo  
caxis([0 0.01])                                   % Eixo do valor de potência (ajustar)
xlabel 'Tempo (s)'
ylabel 'Frequência (Hz)'
colorbar
set(gca, 'Fontsize', 10)

subplot(4,4,[13 14])
plot(F,log(mean(Psd1)),'k','linewidth',1)        % log PSD - AMB
hold on
plot(F,log(mean(Psd2)),'r','linewidth',1)        % log PSD - EXP
xlim([0 40])                                     % Eixo das frequências   
plot(F,log(mean(Psd1)+std(Psd1)/sqrt(5)),'k--')  % log EPM positivo do AMB
plot(F,log(mean(Psd1)-std(Psd1)/sqrt(5)),'k--')  % log EPM negativo do AMB
plot(F,log(mean(Psd2)+std(Psd2)/sqrt(5)),'r--')  % log EPM positivo do EXP
plot(F,log(mean(Psd2)-std(Psd2)/sqrt(5)),'r--')  % log EPM negativo do EXP
title 'AMB vs. EXP'
xlabel 'Frequência (Hz)'
ylabel 'Potência (log)'
Legend={'AMB.','EXP.'};
legend(Legend,'box','off')
box off
set(gca, 'Fontsize', 10)

subplot(4,4,[15 16])
bar([1:1:5],[mean(NormAMB)],0.2,'k')        % Apresenta 5 barras da potência normalizada de AMB
hold on
bar([1.25:1:5.25],[mean(NormEXP)],0.2,'r')  % Apresenta 5 barras da potência normalizada de EXP
errorbar([1:1:5],[mean(NormAMB)],[std(NormAMB)/sqrt(5)],'.k')       %erro padrão da média AMB
errorbar([1.25:1:5.25],[mean(NormEXP)],[std(NormEXP)/sqrt(5)],'.k') %erro padrão da média EXP
xlabel 'Banda da frequência (Hz)'
ylabel 'Potência norm'
ylim([0.6 1.4])
xlim([0.75 5.75])
box off
set(gca, 'Fontsize', 10)
set(gca,'xtick',[1.125:1:5.125],'xticklabels',{'2-4','7-12','15-25','29-33','36-40'}) %opçao de declaração numérica
%set(gca,'xtick',[1:5],'xticklabels',{'Delta','Theta','Beta','Gamma'})                %opção de declaração nominal  
for b=1:5 % número de bandas
    if p(b)<0.05
        text(b,1.3,'*')  % se valor de p for significativo, será reportado
    end
end


%% Coerência  - avaliar coerência entre dois núcleos
% Os sinais já devem estar filtrados e limpos
% Vamos considerar os sinais pronto X_F e Y_F
clear LFP*
LFP1 = X_F; % núcleo 1
LFP2 = Y_F; % núcleo 2

dt = 1/srate;
time_vector = dt:dt:length(LFP1(1,:))/srate;

subplot(311)
plot(time_vector,mean(LFP1,1),'k-')
hold on
plot(time_vector,mean(LFP2,1)-1,'b-')
title(['LFPs ',' Núcleos'],'fontsize',14)
xlim([0 600])
ylim([-4 1])
ylabel('mV')
xlabel('Tempo (s)')
set(gca,'FontSize',12)


subplot(312)
plot(time_vector,eegfilt(mean(LFP1,1),srate,2,6),'k-') % filtra as bandas de acordo com o interesse de observar coerência
hold on
plot(time_vector,eegfilt(mean(LFP2,1),srate,2,6),'b-') % filtra as bandas de acordo com o interesse de observar coerência
xlim([290 310])
ylabel('mV')
title('LFP Filtrado','fontsize',13)
xlabel('Time (s)')
set(gca,'FontSize',12)


% Calculando o Espectro de Coerencias
window = 5*srate;   % tamanho da janela n*srate
overlap = window/2; % sobreposição da janela
%Avalia o quadrado da magnitude de coerência (O valor absoluto ao longo dos pontos de frequência no tempo)
[Cxy, F] = mscohere(mean(LFP2,1),mean(LFP1,1),window,overlap,2^13,srate); 

subplot(313)
plot(F,Cxy,'k')
xlim([0 50])
xlabel('Frequência (Hz)')
ylabel(' Coerência','fontsize',13)
ylim([0 1])
title('Coerência')
set(gcf,'color','white')
set(gca,'FontSize',12)