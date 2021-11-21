%% Setup
num_neurons = 5; % Number of neurons

sim_time = 35;
tstep = 1e-2;
init_conditions = zeros(num_neurons, 1);
refractory_times = zeros(num_neurons, 1);
offset_currents = 1.1*ones(num_neurons, 1);
synaptic_density = zeros(num_neurons, num_neurons)+0.1;
noise_amplitude = randn(num_neurons, 1);

sum_neurons = cell(num_neurons); % zeros(sim_time/tstep, num_neurons);


% Weight matrix variables
mean = [-10 -5 -1 0 1 5 10];
variance = [0 0.1 1 5 10];

% Final data output
data_matrix = zeros(size(mean,2) * size(variance,2), 3); % num_combinations x 3

%% Main Loop
row_num = 1
for i = mean
    for j = variance
        W = i + j * randn(num_neurons);
        W(1:1+size(W,1):end) = 0;
        [spkFull, NetParams, V] = SimLIFNet( ...
            W, ...
            'simTime',sim_time, ...
            'tstep',tstep, ... 
            'offsetCurrents',1.1*(zeros(num_neurons, 1)+1), ...
            'synapticDensity',synaptic_density, ...
            'initialConditions',init_conditions, ... % Membrane voltage
            'refractoryTime',refractory_times, ... 
            'offsetCurrents',offset_currents, ...
            'noiseAmplitude',noise_amplitude, ...
            'displayProgress', 0, ... % Suppress output
            'plotResults', 0 ... % Suppress output
        );
        
        % Add summed voltages to cell array
        % sum_neurons(find(mean==i), find(mean==j)) = {sum(V)};
        data_matrix(row_num, 1) = i; % mean
        data_matrix(row_num, 2) = j;
        data_matrix(row_num, 3) = std(sum(V));
        
        
        
        row_num = row_num + 1;
    end
end

%% Outputting Data Matrix
writematrix(data_matrix, 'data.csv');

%% Plotting
for i = 1:length(mean)
    figure (i) % One window per mean
    hold on
    
    % Overlay all variances on each other
    for j = 1:length(variance)
        plot(sum_neurons{i, j});
    end
    
    hold off
end

%% Test code for running the sim
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