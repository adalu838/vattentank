%% Laboration 3 i TSFS06: Diagnos & Övervakning
clear all; close all; clc;
%% ==================================================
%  Sätt modellparametrar för simulering
%  ==================================================

% Parametrar som används för Simulering av vattentanken
% d1--d6 Parametrar som beskriver dynamiken för vattentankssystemet.
%        Se labbkompendiet för detaljer.
% Uref   Referenssignal till regulator.
% t      Tidsvektor för referenssignal.
% h1Init Initialvärde för vattennivå i tank 1. Default är Uref(1).
% h2Init Initialvärde för vattennivå i tank 2. Beräknas utifrån
%        stationär punkt med avseende på h1Init.
%
load lab3init

h1Init=Uref(1,1);        % Initialnivå i tank 1
h2Init=(d3/d4)^2*h1Init; % Initialnivå i tank 2 
d = [d1 d2 d3 d4 d5 d6]; % Spara alla modellparametrar i en vektor

watertankparams.x0 = [h1Init; h2Init]; % Initialnivån i resp. tank
watertankparams.d  = d; % Modellparametrar

% Sätt slumptalsfrö och brusintensiteter för simulering av mätbrus
NoiseSeed = floor(abs(randn(1,4)*100));
NoiseCov  = 1e-4*[5 0.25 0.5 0.5];

%% ==================================================
%  Design av residualgenerator 'obs1'
%  via observatörsdesign
%  ==================================================

% Linjärisering av vattentankssystemet för h1=4:
Gsys=tanklinj(4, d);
% Anpassa A och C-matrisen så att de passar
A  = Gsys.a(1,1);
C  = Gsys.c(1,1);

P  = [-1]; % Placering av polerna
K1 = obsgain(A,C,P);

% Spara parametrarna som skickas in till observatören i 
obs1params.x0 = h1Init; % Initialvärde på observatörens tillstånd
obs1params.K1 = K1;     % Observatörsförstärkningen
obs1params.d  = d;      % Modellparametrar

%% ==================================================
%  Design av residualgenerator consrel1
%  via konsistensrelation med dynamik
%  ==================================================

% Sätt parametrar för konsistensrelation
c1params.x0   = -h1Init; % Initialvärde för residualgeneratorns tillstånd
c1params.alfa = 2;       % Placering av polen i -alfa
c1params.d    = d;       % Modellparametrar


%% ==================================================
%  Tröskelsättning
%  ==================================================

Jnorm=ones(1,3); % Default är alla trösklar satta till 1

%% ==================================================
%  Simulera systemet
%  simuleringen kan antingen göras genom att välja menyn
%  Simulation->Start i Simulink fönstret 
%  eller exekvera nedanstående rad
%  ==================================================

sim('TSFS06Lab3');

%% ==================================================
%  Definiera beslutsstrukturen via s0 och s1
%  Felfria fallet NF ska stå först
%  ==================================================

% Beslut för residualer under tröskeln
%s0 = ones(3,11);
s0 = ones(3,7); %förenklad variant

% Beslut för residualer över tröskeln
%s1 = zeros(3,11);
s1 = zeros(3,7); %förenklad variant

%% ==================================================
%  Beräkna diagnoser under ett enkelfelsantagande
%  ==================================================
[S,alarm] = decisioncalc(T,s0,s1);


%% ==================================================
%  Plotta resultatet
%  ==================================================

% Förslag på plottar
figure(1)
plot( tut, y )
figure(2)
plot( tut, alarm )
figure(3)
plot( tut, res )
figure(4)
plot( tut,T )

figure(5)
% Kräver att felmoderna är definierade i samma ordning
% i 'S' som i 'name'.
name={'NF', 'Fa', 'Fh2', 'Ff1', 'Fl2', 'Fl3', 'Fc1'};
%name={'NF', 'Fa', 'Fh1', 'Fh2', 'Ff1', 'Ff2', 'Fl1', 'Fl2', 'Fl3',...
%'Fc1', 'Fc2'};

% Plottar diagnosbeslutet för de olika felmoderna enligt S 
% och namnger dem efter name.
for n=1:length(name)
  subplot(3,3,n)
  plot(tut,S(:,n))
  title(name{n})
  axis([min(tut) max(tut) -0.1 1.1])
end
