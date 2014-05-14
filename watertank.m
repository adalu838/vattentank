function watertank(block)

  setup(block);

function setup(block)
%% ===============================================================
%% Definiera modellparametrar ====================================
%% ===============================================================
  numstates = 2; % Antal kontinuerliga tillstånd i funktionen
%% ===============================================================
  
  numparams = 1;

  %% Register number of dialog parameters   
  block.NumDialogPrms = numparams;

  %% Register number of input and output ports
  block.NumInputPorts  = 2;
  block.NumOutputPorts = 4;

  %% Setup functional port properties to dynamically
  %% inherited.
  block.SetPreCompInpPortInfoToDynamic;
  block.SetPreCompOutPortInfoToDynamic;
 
  block.InputPort(1).Dimensions        = 1;
  block.InputPort(2).Dimensions        = 10;
  block.InputPort(1).DirectFeedthrough = true;
  block.OutputPort(1).Dimensions       = 1;
  block.OutputPort(2).Dimensions       = 1;
  block.OutputPort(3).Dimensions       = 1;
  block.OutputPort(4).Dimensions       = 1;
  
  %% Set block sample time to continuous
  block.SampleTimes = [0 0];
  
  %% Setup Dwork
  block.NumContStates = numstates;

  %% Register methods
  block.RegBlockMethod('InitializeConditions',    @InitConditions);  
  block.RegBlockMethod('Outputs',                 @Output);  
  block.RegBlockMethod('Derivatives',             @Derivative);  
  block.RegBlockMethod('SetInputPortSamplingMode', @SetInpPortFrameData);
     
%endfunction

function InitConditions(block)
%% ===============================================================
%% Definiera initialtillståndet ==================================
%% ===============================================================
  x0 = block.DialogPrm(1).Data.x0;
%% ===============================================================

  %% Initialize Dwork
  block.ContStates.Data = x0;
  
%endfunction

%% ===============================================================
%% Definiera utsignalsekvationen =================================
%% ===============================================================
function Output(block)
  d  = block.DialogPrm(1).Data.d;
  h1 = block.ContStates.Data(1);
  h2 = block.ContStates.Data(2);

  fh1 = block.InputPort(2).Data(2);
  fh2 = block.InputPort(2).Data(3);
  ff1 = block.InputPort(2).Data(4);
  ff2 = block.InputPort(2).Data(5);
  fl1 = block.InputPort(2).Data(6);
  fl3 = block.InputPort(2).Data(8);
  fc1 = block.InputPort(2).Data(9);
  fc2 = block.InputPort(2).Data(10);

  y1 = h1 + fh1;
  y2 = h2 + fh2;
  y3 = d(5)*(1-fc1)*(1-fl1)*sqrt(h1) + ff1;
  y4 = d(6)*(1-fc2)*(1-fl3)*sqrt(h2) + ff2;
  
  block.OutputPort(1).Data = y1;
  block.OutputPort(2).Data = y2;
  block.OutputPort(3).Data = y3;
  block.OutputPort(4).Data = y4;
  
%endfunction

%% ===============================================================
%% Definiera de dynamiska ekvationerna ===========================
%% ===============================================================
function Derivative(block)
  d  = block.DialogPrm(1).Data.d;

  h1 = block.ContStates.Data(1);
  h2 = block.ContStates.Data(2);

  u = block.InputPort(1).Data(1);

  fa  = block.InputPort(2).Data(1);
  fl1 = block.InputPort(2).Data(6);
  fl2 = block.InputPort(2).Data(7);
  fc1 = block.InputPort(2).Data(9);
  fc2 = block.InputPort(2).Data(10);
  
  dh1 = d(1)*u-d(2)*(1-fc1)*sqrt(h1)+fa;  
  dh2 = d(3)*(1-fc1)*(1-fl1)*(1-fl2)*sqrt(h1)-d(4)*(1-fc2)*sqrt(h2);

  block.Derivatives.Data(1) = dh1;
  block.Derivatives.Data(2) = dh2;

function SetInpPortFrameData(block, idx, fd)
  
  block.InputPort(idx).SamplingMode = fd;
  block.OutputPort(1).SamplingMode  = fd;
  block.OutputPort(2).SamplingMode  = fd;
  block.OutputPort(3).SamplingMode  = fd;
  block.OutputPort(4).SamplingMode  = fd;
