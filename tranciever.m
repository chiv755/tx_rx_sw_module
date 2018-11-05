%% Main source file to run Simualtion of traceiver
%
%   <put short description/help> here
%
%  Author: Chugunov IV
%  Created: 28.10.2018
%  Modified: 28.10.2018
%
%addpath C:\Users\Ilya\Desktop\inst\taf;


clear all;clc; close all;
%% configuration 

config = init_config();

%% data to tranmit generation
data = randi(2, 1, config.numBits);
data = data - ones(size(data)); 

%% TX (transceiver) simulation
txOut = tx(data, config);

%% Channel ( imperfections) simulation
channelOut =  add_noise(txOut, config.awgn); % add white gaussian noise

%% RX (receiver) simulation
rxOut = rx(channelOut, config);

%% Plot results
plot_results(data, rxOut, config)