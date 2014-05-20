function consrel7(block)

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

  %% Initialize Dwork
  block.ContStates.Data = x0;
  
%endfunction

%% ===============================================================
%  Definiera utsignalsekvationen
%  ===============================================================
function Output(block)
  w  = block.ContStates.Data;
  y2 = block.InputPort(1).Data(3);
  r  = w+y2;
  
  block.OutputPort(1).Data = r;
  
%endfunction

%% ===============================================================
%  Definiera de dynamiska ekvationerna
%  ===============================================================
function Derivative(block)
alfa = block.DialogPrm(1).Data.alfa;
d    = block.DialogPrm(1).Data.d;

y3 = block.InputPort(1).Data(4);
y2 = block.InputPort(1).Data(3);
w  = block.ContStates.Data;

dw = -alfa*(w+y2)-d(3)/d(5)*y3+d(4)*sqrt(max(0,y2));

block.Derivatives.Data = dw;
