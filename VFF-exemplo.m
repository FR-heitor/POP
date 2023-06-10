%% O código seguinte é uma adaptação do código de Gonzales et al. 2020 %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Phase-locking value - slow phase modulation - Manager program
% 
% This program loads and preprocesses all inputs necessary to run 
% the core function plv_modindex.m (Gonzalez et al. 2020, 
% "Communication through coherence by means of cross-frequency coupling"
% https://doi.org/10.1101/2020.03.09.984203)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Carregar banco de dados
clear plv_modindex_comodulogram PLV1

%Condição 1
clear lfp*
for a = % número de animais
lfp1 = X_1(a,t1*srate:t2*srate); % Sinal 1
lfp2 = Y_1(a,t1*srate:t2*srate); % Sinal 2

addpath 'G:\\Phase-Locking-Value---Modulation-Index-master'
% Definindo Parâmetros
srate = 600; % taxa de amostragem
slow_vector = 0:1:10; % vetor de frequências para filtrar as ondas lentas
fast_vector = 10:5:50; % vetor de frequências para filtrar as ondas rápidas
slow_BandWidth = 3; % largura de banda do filtro (ondas lentas)
fast_BandWidth = 10; % largura de banda do filtro (ondas rápidas)
numbin = 20; % número de divisões de fase (ondas lentas)
data_length = length(lfp1); % obtendo o comprimento da série temporal lfp

% Pré-alocação
FastPhase1 = zeros(length(fast_vector), data_length); % pré-alocação
FastPhase2 = zeros(length(fast_vector), data_length); % pré-alocação
SlowPhase = zeros(length(slow_vector), data_length); % pré-alocação

% Obtendo a série temporal da fase de frequência rápida
for ii=1:length(fast_vector)
    Ff1 = fast_vector(ii); % selecionando frequência (corte baixo)
    Ff2=Ff1+fast_BandWidth; % selecionando frequência (corte alto)
    FastFreq1=eegfilt(lfp1,srate,Ff1,Ff2); % filtrando lfp1
    FastFreq2=eegfilt(lfp2,srate,Ff1,Ff2); % filtrando lfp2
    FastPhase1(ii, :) = angle(hilbert(FastFreq1)); % obtendo a fase instantânea de lfp1
    FastPhase2(ii, :) = angle(hilbert(FastFreq2)); % obtendo a fase instantânea de lfp2
end
 
% Obtendo a série temporal da fase de frequência lenta
for jj=1:length(slow_vector)
    Sf1 = slow_vector(jj); % selecionando frequência (corte baixo)
    Sf2 = Sf1 + slow_BandWidth; % selecionando frequência (corte alto)
    SlowFreq = eegfilt(lfp1,srate,Sf1,Sf2); % filtrando lfp com a referência de onda lenta
    SlowPhase(jj, :) = angle(hilbert(SlowFreq)); % obtendo a fase instantânea da onda lenta
end

% Loop através das frequências e calcula o comodulograma plv_modindex
plv_modindex_comodulogram = zeros(size(FastPhase1,1),size(SlowPhase,1)); % pré-alocação
for i = 1:size(SlowPhase,1) % loop através das frequências lentas
    for j = 1:size(FastPhase1,1) % loop através das frequências rápidas   
        plv_phase_modindex = plv_modindex(FastPhase1(j,:)',FastPhase2(j,:)',SlowPhase(i,:)',numbin); % cálculo do plv_modindex
        plv_modindex_comodulogram(j,i) = plv_phase_modindex; % armazenando os resultados na variável plv_modindex_comodulogram
    end  
end
PLV1(:,:,a)  = plv_modindex_comodulogram;
end

%Condição 2
clear lfp*
for a = % número de animais
lfp1 = X_2(a,t1*srate:t2*srate); % Sinal 1
lfp2 = Y_2(a,t1*srate:t2*srate); % Sinal 2

addpath 'G:\\Phase-Locking-Value---Modulation-Index-master'
% Definindo Parâmetros
srate = 600; % taxa de amostragem
slow_vector = 0:1:10; % vetor de frequências para filtrar as ondas lentas
fast_vector = 10:5:50; % vetor de frequências para filtrar as ondas rápidas
slow_BandWidth = 3; % largura de banda do filtro (ondas lentas)
fast_BandWidth = 10; % largura de banda do filtro (ondas rápidas)
numbin = 20; % número de divisões de fase (ondas lentas)
data_length = length(lfp1); % obtendo o comprimento da série temporal lfp

% Pré-alocação
FastPhase1 = zeros(length(fast_vector), data_length); % pré-alocação
FastPhase2 = zeros(length(fast_vector), data_length); % pré-alocação
SlowPhase = zeros(length(slow_vector), data_length); % pré-alocação

% Obtendo a série temporal da fase de frequência rápida
for ii=1:length(fast_vector)
    Ff1 = fast_vector(ii); % selecionando frequência (corte baixo)
    Ff2=Ff1+fast_BandWidth; % selecionando frequência (corte alto)
    FastFreq1=eegfilt(lfp1,srate,Ff1,Ff2); % filtrando lfp1
    FastFreq2=eegfilt(lfp2,srate,Ff1,Ff2); % filtrando lfp2
    FastPhase1(ii, :) = angle(hilbert(FastFreq1)); % obtendo a fase instantânea de lfp1
    FastPhase2(ii, :) = angle(hilbert(FastFreq2)); % obtendo a fase instantânea de lfp2
end
 
% Obtendo a série temporal da fase de frequência lenta
for jj=1:length(slow_vector)
    Sf1 = slow_vector(jj); % selecionando frequência (corte baixo)
    Sf2 = Sf1 + slow_BandWidth; % selecionando frequência (corte alto)
    SlowFreq = eegfilt(lfp1,srate,Sf1,Sf2); % filtrando lfp com a referência de onda lenta
    SlowPhase(jj, :) = angle(hilbert(SlowFreq)); % obtendo a fase instantânea da onda lenta
end

% Loop através das frequências e calcula o comodulograma plv_modindex
plv_modindex_comodulogram = zeros(size(FastPhase1,1),size(SlowPhase,1)); % pré-alocação
for i = 1:size(SlowPhase,1) % loop através das frequências lentas
    for j = 1:size(FastPhase1,1) % loop através das frequências rápidas   
        plv_phase_modindex = plv_modindex(FastPhase1(j,:)',FastPhase2(j,:)',SlowPhase(i,:)',numbin); % cálculo do plv_modindex
        plv_modindex_comodulogram(j,i) = plv_phase_modindex; % armazenando os resultados na variável plv_modindex_comodulogram
    end  
end
PLV2(:,:,a)  = plv_modindex_comodulogram;
end