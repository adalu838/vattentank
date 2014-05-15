function obs4(block)

  setup(block);

function setup(block)
%% ===============================================================
%  Definiera parametrar för residualgeneratorn 
%  ===============================================================
  numstates = 2; % Antal kontinuerliga tillstånd i funktionen
%  ===============================================================

  numparams = 1;

  %% Register number of dialog parameters   
  block.NumDialogPrms = numparams;

  %% Register number of input and output ports
  block.NumInputPorts  = 1;
  block.NumOutputPorts = 1;

  %% Setup functional port properties to dynamically
  %% inherited.
  block.SetPreCompInpPortInfoToDynamic;
  block.SetPreCompOutPortInfoToDynamic;
 
  block.InputPort(1).Dimensions        = 5;
  block.InputPort(1).DirectFeedthrough = true;
  block.OutputPort(1).Dimensions       = 1;
  
  %% Set block sample time to continuous
  block.SampleTimes = [0 0];
  
  %% Setup Dwork
  block.NumContStates = numstates;

  %% Register methods
  block.RegBlockMethod('InitializeConditions',    @InitConditions);  
  block.RegBlockMethod('Outputs',                 @Output);  
  block.RegBlockMethod('Derivatives',             @Derivative);  
  
%endfunction

function InitConditions(block)
%% ===============================================================
%  Definiera initialtillståndet
%  ===============================================================
  x0 = block.DialogPrm(1).Data.x0;
%  ===============================================================
  block.ContStates.Data = x0;
  
%endfunction

%% ===============================================================
%  Definiera utsignalsekvationen
%  ===============================================================
function Output(block)
  x2hat = block.ContStates.Data(2);
  y4    = block.InputPort(1).Data(5);
  d  = block.DialogPrm(1).Data.d;

  r     = y4-d(6)*sqrt(max(0,x2hat));
  block.OutputPort(1).Data = r;
  
%endfunction

%% ===============================================================
%  Definiera de dynamiska ekvationerna
%  ===============================================================
function Derivative(block)
K4 = block.DialogPrm(1).Data.K;
d  = block.DialogPrm(1).Data.d;

u     = block.InputPort(1).Data(1);
y4    = block.InputPort(1).Data(5);
x1hat = block.ContStates.Data(1);
x2hat = block.ContStates.Data(2);

dx1hat = d(1)*u-d(2)*sqrt(max(0,x1hat))+K4(1)*(y4-d(6)*sqrt(max(0,x2hat)));
dx2hat = d(3)*sqrt(max(0,x1hat))-d(4)*sqrt(max(0,x2hat))+K4(2)*(y4-d(6)*sqrt(max(0,x2hat)));

block.Derivatives.Data = [dx1hat dx2hat];
