% Setup
num_neurons = 5; % Number of neurons
sim_time = 35;
tstep = 1e-2;
v_neurons = zeros(sim_time/tstep, num_neurons);

% Weight matrix variables
mean = 5;%[-10 -5 -1 0 1 5 10];
variance = 2;%[-10 -5 -1 0 1 5 10];

for i = mean
    for j = variance
        W = i + j * randn(num_neurons);
        W(1:1+size(W,1):end) = 0;
        [spkFull, NetParams, V] = SimLIFNet(W,'simTime',sim_time,'tstep',tstep,... 
          'offsetCurrents',1.1*(zeros(N, 1)+1),'synapticDensity',zeros(N, N)+0.1,...
          'initialConditions',zeros(N, 1));
    end
end
% surf(W)

% Test code for running the sim
% W = [0 -0.2; -0.2 0];
% % Asynchronous behaviour at low driving currents (1.1) and synaptic
% % densities of 3, and initial membrane potentials of 0.4 and 0
% [spkFull NetParams V] = SimLIFNet(W,'simTime',35,'tstep',1e-2,... 
%           'offsetCurrents',1.1*[1; 1],'synapticDensity',[3; 3],...
%           'initialConditions',[0.4; 0]);
% % Now observe synchronous behavior with stronger forcing, 1.6
%  [spkFull NetParams V] = SimLIFNet(W,'simTime',35,'tstep',1e-2,... 
%           'offsetCurrents',1.6*[1; 1],'synapticDensity',[3; 3],...
%           'initialConditions',[0.4; 0]);