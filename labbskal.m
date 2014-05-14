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
%  Design av residualgenerator 'obs1'
%  via observat�rsdesign
%  ==================================================

% Linj�risering av vattentankssystemet f�r h1=4:
Gsys=tanklinj(4, d);
% Anpassa A och C-matrisen s� att de passar
A  = Gsys.a(1,1);
C  = Gsys.c(1,1);

P  = [-1]; % Placering av polerna
K1 = obsgain(A,C,P);

% Spara parametrarna som skickas in till observat�ren i 
obs1params.x0 = h1Init; % Initialv�rde p� observat�rens tillst�nd
obs1params.K1 = K1;     % Observat�rsf�rst�rkningen
obs1params.d  = d;      % Modellparametrar

%% ==================================================
%  Design av residualgenerator consrel1
%  via konsistensrelation med dynamik
%  ==================================================

% S�tt parametrar f�r konsistensrelation
c1params.x0   = -h1Init; % Initialv�rde f�r residualgeneratorns tillst�nd
c1params.alfa = 2;       % Placering av polen i -alfa
c1params.d    = d;       % Modellparametrar


%% ==================================================
%  Tr�skels�ttning
%  ==================================================

Jnorm=ones(1,3); % Default �r alla tr�sklar satta till 1

%% ==================================================
%  Simulera systemet
%  simuleringen kan antingen g�ras genom att v�lja menyn
%  Simulation->Start i Simulink f�nstret 
%  eller exekvera nedanst�ende rad
%  ==================================================

sim('TSFS06Lab3');

%% ==================================================
%  Definiera beslutsstrukturen via s0 och s1
%  Felfria fallet NF ska st� f�rst
%  ==================================================

% Beslut f�r residualer under tr�skeln
%s0 = ones(3,11);
s0 = ones(3,7); %f�renklad variant

% Beslut f�r residualer �ver tr�skeln
%s1 = zeros(3,11);
s1 = zeros(3,7); %f�renklad variant

%% ==================================================
%  Ber�kna diagnoser under ett enkelfelsantagande
%  ==================================================
[S,alarm] = decisioncalc(T,s0,s1);


%% ==================================================
%  Plotta resultatet
%  ==================================================

% F�rslag p� plottar
figure(1)
plot( tut, y )
figure(2)
plot( tut, alarm )
figure(3)
plot( tut, res )
figure(4)
plot( tut,T )

figure(5)
% Kr�ver att felmoderna �r definierade i samma ordning
% i 'S' som i 'name'.
name={'NF', 'Fa', 'Fh2', 'Ff1', 'Fl2', 'Fl3', 'Fc1'};
%name={'NF', 'Fa', 'Fh1', 'Fh2', 'Ff1', 'Ff2', 'Fl1', 'Fl2', 'Fl3',...
%'Fc1', 'Fc2'};

% Plottar diagnosbeslutet f�r de olika felmoderna enligt S 
% och namnger dem efter name.
for n=1:length(name)
  subplot(3,3,n)
  plot(tut,S(:,n))
  title(name{n})
  axis([min(tut) max(tut) -0.1 1.1])
end
