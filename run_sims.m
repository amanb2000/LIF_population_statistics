%% Setup
num_neurons = 5; % Number of neurons

sim_time = 35;
tstep = 1e-2;
init_conditions = zeros(num_neurons, 1);
refractory_times = zeros(num_neurons, 1);
offset_currents = 1.1*ones(num_neurons, 1);
synaptic_density = zeros(num_neurons, num_neurons)+0.1;
noise_amplitude = .01*randn(num_neurons, 1);

sum_neurons = cell(num_neurons); % zeros(sim_time/tstep, num_neurons);

num_repeats = 15;

% Weight matrix variables
% mean = [-10 -5 -1 0 1 5 10];
% variance = [0 0.1 1 5 10];
n_mean = 30;
n_var = 10;

var_max = 1
mean_range = 1

means = 1:n_mean
variance = 1:n_var
means = (means*(1/n_mean))*2*mean_range - mean_range
variance = (variance*(1/n_var))*var_max

% Final data output
num_experiments = size(means,2)*size(variance,2)*num_repeats;
data_matrix = zeros(size(means,2)*size(variance,2), 5); % num_combinations x 4

%% Main Loop
row_num = 1
exp_num = 1
for i = means
    for j = variance
        mean_max_spike=0;
        mean_synch = 0;
        
        synchs = zeros(1,num_repeats);
        
        for k = 1:num_repeats
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



            % Calculating how many times the maximally active neuron fired.
            max_spike = -1;
            for nnum = 1:num_neurons
                if size(spkFull{nnum},2) > max_spike
                    max_spike = size(spkFull{nnum},2);
                end
            end

            mean_max_spike = mean_max_spike + max_spike/num_repeats;
%             mean_synch = mean_synch + std(sum(V))/num_repeats;
               
            synchs(k) = std(sum(V))/num_repeats;
            
            fprintf("%d of %d\n",exp_num, num_experiments);
            exp_num = exp_num + 1;
        end
        
        % Add summed voltages to cell array
        % sum_neurons(find(mean==i), find(mean==j)) = {sum(V)};
        data_matrix(row_num, 1) = i; % mean
        data_matrix(row_num, 2) = j;
        data_matrix(row_num, 3) = mean_max_spike;
        data_matrix(row_num, 4) = mean(synchs);
        data_matrix(row_num, 5) = std(synchs);
        
        row_num = row_num + 1;
    end
end

%% Outputting Data Matrix
writematrix(data_matrix, 'data_2.csv');

%% Test code for arbitrary condition

arb_mean = 0.01
arb_var = 0.01

W = arb_mean + arb_var * randn(num_neurons);
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
    'displayProgress', 1, ... % Suppress output
    'plotResults', 1 ... % Suppress output
);

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