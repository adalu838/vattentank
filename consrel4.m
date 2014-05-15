function consrel4(block)

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
  y1 = block.InputPort(1).Data(2);
  r  = w+y1;
  
  block.OutputPort(1).Data = r;
  
%endfunction

%% ===============================================================
%  Definiera de dynamiska ekvationerna
%  ===============================================================
function Derivative(block)
alfa = block.DialogPrm(1).Data.alfa;
d    = block.DialogPrm(1).Data.d;

u  = block.InputPort(1).Data(1);
y1 = block.InputPort(1).Data(2);
w  = block.ContStates.Data;

dw = -alfa*(w+y1)-d(1)*u+d(2)*sqrt(max(0,y1));

block.Derivatives.Data = dw;
