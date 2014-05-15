function obs2(block)

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
  y2    = block.InputPort(1).Data(3);

  r     = y2-x2hat;
  block.OutputPort(1).Data = r;
  
%endfunction

%% ===============================================================
%  Definiera de dynamiska ekvationerna
%  ===============================================================
function Derivative(block)
K2 = block.DialogPrm(1).Data.K;
d  = block.DialogPrm(1).Data.d;

u     = block.InputPort(1).Data(1);
y2    = block.InputPort(1).Data(3);
x1hat = block.ContStates.Data(1);
x2hat = block.ContStates.Data(2);

dx1hat = d(1)*u-d(2)*sqrt(x1hat)+K2(1)*(y2-x2hat);
dx2hat = d(3)*sqrt(x1hat)-d(4)*sqrt(x2hat)+K2(2)*(y2-x2hat);

block.Derivatives.Data = [dx1hat dx2hat];
