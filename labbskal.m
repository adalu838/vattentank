%% Laboration 3 i TSFS06: Diagnos & �vervakning
clear all; close all; clc;
%% ==================================================
%  S�tt modellparametrar f�r simulering
%  ==================================================

% Parametrar som anv�nds f�r Simulering av vattentanken
% d1--d6 Parametrar som beskriver dynamiken f�r vattentankssystemet.
%        Se labbkompendiet f�r detaljer.
% Uref   Referenssignal till regulator.
% t      Tidsvektor f�r referenssignal.
% h1Init Initialv�rde f�r vattenniv� i tank 1. Default �r Uref(1).
% h2Init Initialv�rde f�r vattenniv� i tank 2. Ber�knas utifr�n
%        station�r punkt med avseende p� h1Init.
%
load lab3init

h1Init=Uref(1,1);        % Initialniv� i tank 1
h2Init=(d3/d4)^2*h1Init; % Initialniv� i tank 2 
d = [d1 d2 d3 d4 d5 d6]; % Spara alla modellparametrar i en vektor

watertankparams.x0 = [h1Init; h2Init]; % Initialniv�n i resp. tank
watertankparams.d  = d; % Modellparametrar

% S�tt slumptalsfr� och brusintensiteter f�r simulering av m�tbrus
NoiseSeed = floor(abs(randn(1,4)*100));
NoiseCov  = 1e-4*[5 0.25 0.5 0.5];

%% ==================================================
%  Design av residualgenerator
%  via observat�rsdesign
%  ==================================================

% Linj�risering av vattentankssystemet f�r h1=4:
Gsys=tanklinj(4, d);

% obs1
% Anpassa A och C-matrisen s� att de passar
A1  = Gsys.a(1,1);
C1  = Gsys.c(1,1);

P1  = [-0.5]; % Placering av polerna
K1 = obsgain(A1,C1,P1);

% Spara parametrarna som skickas in till observat�ren  
obs1params.x0 = h1Init; % Initialv�rde p� observat�rens tillst�nd
obs1params.K = K1;     % Observat�rsf�rst�rkningen
obs1params.d  = d;      % Modellparametrar

% obs2
% Anpassa A och C-matrisen s� att de passar
A2  = Gsys.a(1:2,1:2);
C2  = Gsys.c(2,:);

P2  = [-1 -2]; % Placering av polerna
K2 = obsgain(A2,C2,P2);

% Spara parametrarna som skickas in till observat�ren  
obs2params.x0 = [h1Init h2Init]; % Initialv�rde p� observat�rens tillst�nd
obs2params.K = K2;     % Observat�rsf�rst�rkningen
obs2params.d  = d;      % Modellparametrar

% obs3
% Anpassa A och C-matrisen s� att de passar
A3  = Gsys.a(1,1);
C3  = Gsys.c(3,1);

P3  = [-0.5]; % Placering av polerna
K3 = obsgain(A3,C3,P3);

% Spara parametrarna som skickas in till observat�ren  
obs3params.x0 = h1Init; % Initialv�rde p� observat�rens tillst�nd
obs3params.K = K3;     % Observat�rsf�rst�rkningen
obs3params.d  = d;      % Modellparametrar

% obs4
% Anpassa A och C-matrisen s� att de passar
A4  = Gsys.a(1:2,1:2);
C4  = Gsys.c(4,:);

P4  = [-0.15 -0.2]; % Placering av polerna
K4 = obsgain(A4,C4,P4);

% Spara parametrarna som skickas in till observat�ren  
obs4params.x0 = [h1Init h2Init]; % Initialv�rde p� observat�rens tillst�nd
obs4params.K = K4;     % Observat�rsf�rst�rkningen
obs4params.d  = d;      % Modellparametrar

%% ==================================================
%  Design av residualgenerator consrel1
%  via konsistensrelation med dynamik
%  ==================================================

% S�tt parametrar f�r konsistensrelation
c1params.x0   = -h1Init; % Initialv�rde f�r residualgeneratorns tillst�nd
c1params.alfa = 0.2;       % Placering av polen i -alfa
c1params.d    = d;       % Modellparametrar

% S�tt parametrar f�r konsistensrelation
c2params.x0   = -h2Init; % Initialv�rde f�r residualgeneratorns tillst�nd
c2params.alfa = 0.2;       % Placering av polen i -alfa
c2params.d    = d;       % Modellparametrar

% S�tt parametrar f�r konsistensrelation
c3params.x0   = -h1Init; % Initialv�rde f�r residualgeneratorns tillst�nd
c3params.alfa = 0.4;       % Placering av polen i -alfa
c3params.d    = d;       % Modellparametrar

% S�tt parametrar f�r konsistensrelation
c4params.x0   = -h2Init; % Initialv�rde f�r residualgeneratorns tillst�nd
c4params.alfa = 0.8;       % Placering av polen i -alfa
c4params.d    = d;       % Modellparametrar

% S�tt parametrar f�r konsistensrelation
c5params.x0   = -h1Init; % Initialv�rde f�r residualgeneratorns tillst�nd
c5params.alfa = 0.1;       % Placering av polen i -alfa
c5params.d    = d;       % Modellparametrar

% S�tt parametrar f�r konsistensrelation
c6params.x0   = -h2Init; % Initialv�rde f�r residualgeneratorns tillst�nd
c6params.alfa = 0.2;       % Placering av polen i -alfa
c6params.d    = d;       % Modellparametrar

% S�tt parametrar f�r konsistensrelation
c7params.x0   = -h2Init; % Initialv�rde f�r residualgeneratorns tillst�nd
c7params.alfa = 0.2;       % Placering av polen i -alfa
c7params.d    = d;       % Modellparametrar

%% ==================================================
%  Tr�skels�ttning
%  ==================================================

pfa = 0.001;

if  exist('res','var')
    my = mean(res);
    std_res = std(res);
    Jnorm=norminv(1-pfa/2,my,std_res); % Default �r alla tr�sklar satta till 1
else
    faults = [0 0 0 0 0 0 0 0 0 0];
    Jnorm=ones(1,13);
    sim('TSFS06Lab3');
end

% Jnorm(11) = Jnorm(11)*.85;
% Jnorm(7) = Jnorm(7)*.7;
% Jnorm(5) = Jnorm(5)*1.3;

%% ==================================================
%  Simulera systemet
%  simuleringen kan antingen g�ras genom att v�lja menyn
%  Simulation->Start i Simulink f�nstret 
%  eller exekvera nedanst�ende rad
%  ==================================================
% ==================================================
%  Definiera beslutsstrukturen via s0 och s1
%  Felfria fallet NF ska st� f�rst
%  ==================================================
n = size(res,2);

% Beslut f�r residualer under tr�skeln
s0 = ones(n,11);
%s0 = ones(3,7); %f�renklad variant

% Beslut f�r residualer �ver tr�skeln
s1 = zeros(n,11);
s1(1,[2 3 10]) = 1;
s1(2,[2 4 7 8 10 11]) = 1;
s1(3,[2 5 7 10]) = 1;
s1(4,[2 6 7 8 9 10 11]) = 1;
s1(5,[2 3 10]) = 1;
s1(6,[3 4 7 8 10 11]) = 1;
s1(7,[2 3 10]) = 1;
s1(8,[3 4 7 8 10 11]) = 1;
s1(9,[3 5 7 10]) = 1;
s1(10,[4 6 9 11]) = 1;
s1(11,[2 3 5 7]) = 1;
s1(12,[4 7 6 8 9 10]) = 1;
s1(13,[4 5 8 11]) = 1;

hej = [1 2 4 5 8 9 10];

s1ny = s1(:,hej);

%s1 = zeros(3,7); %f�renklad variant

%% ==================================================
%  Ber�kna diagnoser under ett enkelfelsantagande
%  ==================================================
%         fa' fh1 fh2' ff1' ff2 fl1 fl2' fl3' fc1' fc2
faults = [ 0   0   0    0   0    0   0   0  0   0 ];
sim('TSFS06Lab3');
[S,alarm] = decisioncalc(T,s0,s1);

% ==================================================
%  Plotta resultatet
% ==================================================

% F�rslag p� plottar
figure(1)
plot( tut, y )
figure(2)
plot( tut, alarm )
figure(3)
plot( tut, res )
legend('Obs1','Obs2','Obs3','Obs4','Func1','Func2','Constrel1', ...
       'Constrel2','Constrel3','Constrel4','Constrel5','Constrel6','consrel7');
figure(4)
plot( tut,T )

figure(5)
% Kr�ver att felmoderna �r definierade i samma ordning
% i 'S' som i 'name'.
%name={'NF', 'Fa', 'Fh2', 'Ff1', 'Fl2', 'Fl3', 'Fc1'};
name={'NF', 'Fa', 'Fh1', 'Fh2', 'Ff1', 'Ff2', 'Fl1', 'Fl2', 'Fl3',...
'Fc1', 'Fc2'};

% Plottar diagnosbeslutet f�r de olika felmoderna enligt S 
% och namnger dem efter name.
for n=[1 2 4 5 8 9 10];%1:length(name)
  subplot(3,4,n)
  plot(tut,S(:,n))
  title(name{n})
  axis([min(tut) max(tut) -0.1 1.1])
end
