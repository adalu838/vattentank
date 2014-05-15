function obs3(block)

  setup(block);

function setup(block)
%% ===============================================================
%  Definiera parametrar för residualgeneratorn 
%  ===============================================================
  numstates = 1; % Antal kontinuerliga tillstånd i funktionen
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
  x1hat = block.ContStates.Data;
  y3    = block.InputPort(1).Data(4);
  d     = block.DialogPrm(1).Data.d;

  r     = y3-d(5)*sqrt(x1hat);
  block.OutputPort(1).Data = r;
  
%endfunction

%% ===============================================================
%  Definiera de dynamiska ekvationerna
%  ===============================================================
function Derivative(block)
K3 = block.DialogPrm(1).Data.K;
d  = block.DialogPrm(1).Data.d;

u     = block.InputPort(1).Data(1);
y3    = block.InputPort(1).Data(4);
x1hat = block.ContStates.Data;

dx1hat = d(1)*u-d(2)*sqrt(x1hat)+K3*(y3-d(5)*sqrt(x1hat));

block.Derivatives.Data = dx1hat;
