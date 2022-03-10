fprintf('\n\n\nStarting Steady State PF Run\n\n\n')

diary on
diary(strcat('run_ss',date,'_',num2str(hour(datetime('now'))),'.txt'))
steadyStateCapital
diary off

steadyState_pf_cap

%%


fprintf('\n\n\nStarting No - Trade PF Run\n\n\n')

diary on
diary(strcat('run_cf',date,'_',num2str(hour(datetime('now'))),'.txt'))
cfTransitionCapital
diary off

counterfactualTransition_pf_ca



%% 
fprintf('\n\n\nStarting Trade PF Run\n\n\n')

diary on
diary(strcat('run_ac',date,'_',num2str(hour(datetime('now'))),'.txt'))
actualTransitionCapital
diary off

actualTransition_pf_cap



%% 

fprintf('\n\n\nStarting Offshoring PF Run\n\n\n')

diary on
diary(strcat('run_off',date,'_',num2str(hour(datetime('now'))),'.txt'))
offshoringTransition
diary off

offshoringTransition_pf_ca

%% 

fprintf('\n\n\nStarting Timing PF Run\n\n\n')

diary on
diary(strcat('run_timing',date,'_',num2str(hour(datetime('now'))),'.txt'))
timingTransition
diary off

timingTransition_pf_cap
