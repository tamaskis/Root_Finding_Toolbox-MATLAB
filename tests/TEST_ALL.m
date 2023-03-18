%% TEST_ALL.m
% Root-Finding Toolbox
%
% Runs all unit testing scripts.
%
% Copyright © 2021 Tamas Kis
% Last Update: 2023-03-03
% Website: https://tamaskis.github.io
% Contact: tamas.a.kis@outlook.com
%
% DEPENDENCIES:
%   • Simple Unit Testing Toolbox (https://tamaskis.github.io/Simple_Unit_Testing_Toolbox-MATLAB)
%   • Root-Finding Toolbox (https://tamaskis.github.io/Root_Finding_Toolbox-MATLAB/)



%% SCRIPT SETUP

% clears Workspace and Command Window, closes all figures
clear; clc; close all;

% adds path to testing scripts
addpath('testsuites')

% adds path to toolbox folder
addpath(genpath('../toolbox'))



%% RUNS ALL TESTS

TEST_convergence;
TEST_utils;
TEST_root_bisection;