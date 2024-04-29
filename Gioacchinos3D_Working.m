%% Gioacchino task
% Description of the task to be modelled. Nouns of graspable objects are
% used as stimuli. Participants had to decide if each noun referred to a
% natural or artificat, by performing either a precision or a power reach
% to grasp movement. Response grasp could be compatible or incompatible
% with the grasp typically used to manipulate the objects. They were
% instructed about the grasp at the beginning of the trial for example, preicison with natural objects and
% power with artificial, this didn't change within trial and 50% prob of one
% vs the other.

% There is a trede off between accuracy and velocity, so that one participant
% can choose to be accurate but slower in incongruent trials or quick but 
% inaccurate in incongruent trials 

% Insert habit in the sense that if I see a Big-N I have the habit to grasp
% with a Pw but i might have to overwrite that if I am in the other context
%

% Moreover, the words can suggest a power or precision
% grip (olive vs hammer), hence one can find oneself in a Compatible
% situation (Artifact-Power + Word-power) or in an incompatible
% (Artifact-power + Word-precision) and the same for natural stimuli.

% Instruction about whether artifact or natural should call for a precision
% vs power grip do not change in-trial, the numbers of trial is 512, with 8
% nouns natural + precision, 8 natural + power, 8 artifact + precision, 8
% artifact + power each repeated 16 times

% from now on Pw = power, Pr = precision, A = artifact, N = natural

clear
clear all 
rng('shuffle') 
%% 
% Simulation options after model building below:

% If Sim = 1, simulate single trial. This will reproduce fig. 8. (Although
            % note that, for this and the following simulations, results 
            % will vary each time due to random sampling)

% If Sim = 2, simulate multiple trials where the left context is active 
            % (D{1} = [1 0]'). This will reproduce fig. 10.
             
% If Sim = 3, simulate reversal learning, where the left context is active 
            % (D{1} = [1 0]') in early trials and then reverses in later 
            % trials (D{1} = [0 1]'). This will reproduce fig. 11.
            
% If Sim = 4, run parameter estimation on simulated data with reversal
            % learning. This will reproduce the top panel of fig. 17.
            
% If Sim = 5, run parameter estimation on simulated data with reversal
            % learning from multiple participants under different models
            % (i.e., different parameter values) and perform model comparison. 
            % This will reproduce the bottom panel of fig. 17. This option
            % will also save two structures that include results of model
            % comparison, model fitting, parameter recoverability analyses,
            % and inputs needed for group (PEB) analyses.

Sim = 1;

%% 
%==========================================================================

% D prior beliefs about intial states 
%--------------------------------------------------------------------------

% State that is told to the participant at the beginning of the trial. 

D{1} = [1 0]'; % Context state {'Pw-A + Pr-N','Pw-N + Pr-A'} 


% Behavioural states

D{2} = [1 0 0 0]'; % {'Start', 'Stim', 'Pw', 'Pr'}

% Distriubution for how probable is for a participant to see some stimuli,
% basically saying that is 50% prob of 

D{3} = [.25 .25 .25 .25]'; % Stimuli state {'A-Pw', 'N-Pr', 'N-Pw', 'A-Pr'} 


% A probabilistic likelihood mapping
%--------------------------------------------------------------------------


%---- A{1} STIMULI OBSERVATION
% I start by specifying that all behavioural state lead to a null

for i = 1:4 
    for j = 1:4

    A{1}(:,:,i,j) = [1 1; % Null
                   0 0; % A-Pw
                   0 0; % N-Pw
                   0 0; % A-Pr
                   0 0]; % N-Pr
    end
end 


% Instead Stim behavioral state leads to the observation of one of the four
% stimuli with equal probability in both contexts. However, the probability
% also changes depending on whether it is in Natural or Artificial.
% Basically, when on stimuli (i = 2) and on Natural (j = 1) something
% happens instead when i = 2 and j = 2 something elses happens. 

% or is it 100% to see an A-Pw when you are in the A-pw situation?


%{'Pw-A + Pr-N','Pw-N + Pr-A'} 

Pstim = 0.50;

for i = 2 
    for j = 1

    A{1}(:,:,i,j) = [0 0; % Null
                   Pstim  1-Pstim ; % A-Pw
                   0 0 ; % N-Pw
                   0 0 ; % A-Pr
                   0 0]; % N-Pr
    end
end 

for i = 2 
    for j = 2

    A{1}(:,:,i,j) = [0 0; % Null
                   0 0; % A-Pw
                   1-Pstim  Pstim ; % N-Pw % this is Pw-N
                   0 0 ; % A-Pr
                   0 0 ]; % N-Pr
    end
end 

for i = 2 
    for j = 3

    A{1}(:,:,i,j) = [0 0; % Null
                   0 0 ; % A-Pw
                   0 0 ; % N-Pw
                   1-Pstim  Pstim ; % A-Pr
                   0 0 ]; % N-Pr
    end
end 

for i = 2 
    for j = 4

    A{1}(:,:,i,j) = [0 0; % Null
                   0 0 ; % A-Pw
                   0 0 ; % N-Pw
                   0 0; % A-Pr
                   Pstim  1-Pstim ]; % N-Pr
    end
end 

%---- Now let's map the correct/incorrect observations.

% In all contexts for the start and stim behavioural state the observation
% is null

for i = 1:4
    for j = 1:4
        
        A{2}(:,:,i, j) = [1 1;  % Null 
                          0 0;  % Correct
                          0 0]; % Incorrect
    end
end

% Now I map all the context and grasps to the stimuli 
% {'A-Pw', 'N-Pr', 'N-Pw', 'A-Pr'}


%if I see a A-Pw (j = 1) and I am performing a Pw grasp (i = 3) 
% I am correct in the A-Pw + Pr-N context

for i = 3
    for j = 1
        A{2}(:,:,i, j) = [0 0;  % Null 
                          1 0;  % Correct
                          0 1]; % Incorrect
    end
end

% if I see a N-Pr (j = 2) and I am performing a Pw grasp (i = 3) 

for i = 3
    for j = 2
        A{2}(:,:,i, j) = [0 0;  % Null 
                          0 1;  % Correct
                          1 0]; % Incorrect
    end
end

% if I see a N-Pw (j = 3) and I am performing a Pw grasp (i = 3) I am
% correct in the second context

for i = 3
    for j = 3
        A{2}(:,:,i, j) = [0 0;  % Null 
                          0 1;  % Correct
                          1 0]; % Incorrect
    end
end

% if I see a A-Pr (j = 2) and I am performing a Pw grasp (i = 3) I am always
% incorrect 


for i = 3
    for j = 4
        A{2}(:,:,i, j) = [0 0;  % Null 
                          1 0;  % Correct
                          0 1]; % Incorrect
    end
end

% Let's do the same when performing a PR grasp
% {'A-Pw', 'N-Pr', 'N-Pw', 'A-Pr'}

% If I am seeing a A-Pw (j = 1) and performing a Pr grasp (i = 4) 

for i = 4 
    for j = 1
        A{2}(:,:,i,j) = [0 0;  % Null 
                         0 1;  % Correct
                         1 0]; % Incorrect
   end
end

% If I am seeing a N-Pr (j = 2) and performing a Pr grasp (i = 4) I am
% correct in the first context 

for i = 4 
    for j = 2
        A{2}(:,:,i,j) = [0 0;  % Null 
                         1 0;  % Correct
                         0 1]; % Incorrect
   end
end


% If I am seeing a N-Pw (j = 3) and performing a Pr grasp (i = 4) I am
% correct when in the first context, bc I am following the instruction of
% performing a precision with natural, but incorrect in the second as the
% stimuli would fit. 

for i = 4 
    for j = 3
        A{2}(:,:,i,j) = [0 0;  % Null 
                         1 0;  % Correct
                         0 1]; % Incorrect
   end
end

% If I am seeing a A-Pr (j = 4) and performing a Pr grasp (i = 4) I am
% correct in the second context 

for i = 4 
    for j = 4
        A{2}(:,:,i,j) = [0 0;  % Null 
                         0 1;  % Correct
                         1 0]; % Incorrect
   end
end


% The agent knows its own actions 

%{
for i = 1
    for j = 1:2 
A{3}(:,:,i,j) = [ 1 1;  % Start 
                  0 0;  % Stim
                  0 0;  % Pw
                  0 0]; % Pr
    end 
end

for i = 2
    for j = 1:2 
A{3}(:,:,i,j) = [ 0 0;  % Start 
                  1 1;  % Stim
                  0 0;  % Pw
                  0 0]; % Pr
    end 
end

for i = 3
    for j = 1:2 
A{3}(:,:,i,j) = [ 0 0;  % Start 
                  0 0;  % Stim
                  1 1;  % Pw
                  0 0]; % Pr
    end 
end

for i = 4
    for j = 1:2 
A{3}(:,:,i,j) = [ 0 0;  % Start 
                  0 0;  % Stim
                  0 0;  % Pw
                  1 1]; % Pr
    end 
end

%}

for i = 1:4 
    for j = 1:4 
        A{3}(i,:,i,j) = [1 1];
    end
end



%--------------------------------------------------------------------------
% B Transition probabilities
%--------------------------------------------------------------------------
% Here I need one B for each D I have, therefore, three Bs. Here columns are 
% time t and rows are time t+1. 

% Context is stable within a trial

B{1}(:,:,1) = [1 0;
               0 1];

% Action possibilities  {'Start', 'Stim', 'Pw', 'Pr} this represents the
% only possible set of actions. From start to stim from stim to Pw o Pr e
% from Pw o Pr to start

B{2}(:,:,1) = [0 0 0 0; % Start
               1 0 0 0; % Stim
               0 1 1 1; % Pw
               0 0 0 0]; % Pr

B{2}(:,:,2) = [0 0 0 0;
               1 0 0 0;
               0 0 0 0;
               0 1 1 1];

 
%This is a mtrix to stay in the start state 

B{2}(:,:,3) = [1 1 1 1;
               0 0 0 0;
               0 0 0 0;
               0 0 0 0];


% Stimuli is stable within a trial 

B{3}(:,:,1) = [1 0 0 0;
               0 1 0 0;
               0 0 1 0;
               0 0 0 1];

%--------------------------------------------------------------------------
% C preferences
%--------------------------------------------------------------------------
% Here I specify preferences over each set of outcomes (C), with one matrix
% per outcome modality (numbers of A's). Here columns indicate time points, and rows
% indicate observations (same order as in the corresponding A matrices).
% Therefore columns ALWAYS = T

% The agent starts in start state, sees a stimuli and performs a grasp

T = 3;

No = [size(A{1}, 1) size(A{2},1) size(A{3}, 1)];
% Size gets the first dimension (Row) of A, hence the No matrix is a 5 x 3

% I have 3 A matrices so I will define 3 Cs.

C{1} = zeros(No(1), T); % Stimuli
C{2} = zeros(No(2), T); % Correct incorrect
C{3} = zeros(No(3), T); % Observed behaviour state factor

% Each C maps to its corrsponding A. The agent has no preference about
% being in the start state, seeing the stimuli, performing the grasp. 
% In the second case the agent has prefrence about his grasp, he prefers to
% perform the right grasps depending on the context. He prefers observing a
% grasp that is coherent with the context which I set as Cp (context
% preference)


C{2}(:,:) = [0 0  0;  % Null
             0 0 -5;  % incorrect
             0 0  1];  % correct
                         



% Allowable policies: U or V. 
%==========================================================================

T ;  % number of timesteps defined above in C
Pi = 2; % number of policies
Nf = 3; % number of factors
% U = ones(1,Pi,Nf); % Context and stimuli are not controllable
% V = ones(T-1, Pi, Nf)

V = ones(T-1,Pi, Nf);

% 1. From the start state observe a stimulus
% 2. From the stimulus perform a Pw
% 3. From the stimulus perform a Pr
% 4. From grasp go back to start 

% U(:,:,1) is the context state hence is not controllab.e U(:,:,3) is the 
% stimulus state also not controllable. It can only control the grasp state 
% U(:,:,2).
% {'Start', 'Stim', 'Pw', 'Pr}

V(:,:,2) = [1 1;
            2 2];
            
            
           

% MDP structure 
%==========================================================================

mdp.T = T; % Number of time steps

mdp.V = V; % Allowable policies

mdp.A = A; % state-outcome mapping 

mdp.B = B; % Transition probabilities

mdp.C = C; % preferred states 

mdp.D = D; % prior over initial states



%mdp.a = a; % enable learning over likelihoods probabilities 

%mdp.b = b; % enable learning over state transitions 

%mdp.d = d; % enable learning priors over intial states


mdp.eta = 1; % Learning rate (0-1)

mdp.omega = 1; % Forgettin rate no forgetting (0-1)

mdp.beta = 1; % Lower values increase the influence of habits (E).

mdp.alpha = 8; % Action precision controls the level of randomness.

%--------------------------------------------------------------------------
% Labelling and plotting

% Adding some lables for subsequent plotting 

% One label factor for each D 

label.factor{1} = 'Context State'; label.name{1} = {'Pw-A + Pr-N','Pw-N + Pr-A'};
label.factor{2} = 'Actions'; label.name{2} = {'Start', 'Stim', 'Pw', 'Pr'};
label.factor{3} = 'Stimuli'; label.name{3} = {'A-Pw', 'N-Pr', 'N-Pw', 'A-Pr'};

% One label modeality for each A \
label.modality{1} = 'Condition'; label.outcome{1} = {'Null', 'A-Pw', 'N-Pw', 'A-pr', 'N-Pr'};
label.modality{2} = 'Feedback'; label.outcome{2} = {'Null', 'Correct', 'Incorrect'};
label.modality{3} = 'Observed Behaviour'; label.outcome{3} = {'Start', 'Stim', 'Pw', 'Pr'};


% One label action for each B, usualyl not the first when it says that the
% context does not change. 

label.action{2} =  {'Start', 'Stim', 'Pw', 'Pr'};
mdp.label = label;

% Check whether all matrix dimensions are correct 

mdp = spm_MDP_check(mdp);

if Sim == 1
%% Single trial simulation

% Run active inference estimantion 

MDP = spm_MDP_VB_X_tutorial(mdp);

% We can then use standard plotting routines to visualize simulated neural 
% responses

spm_figure('GetWin','Figure 1'); clf    % display behavior
spm_MDP_VB_LFP(MDP); 

%  and to show posterior beliefs and behavior:

spm_figure('GetWin','Figure 2'); clf    % display behavior
spm_MDP_VB_trial(MDP); 


elseif Sim == 2

    N = 32; %Number of trials 

    MDP = mdp;

    [MDP(1:N)] = deal(MDP);

    MDP = spm_MDP_VB_X_tutorial(MDP);

    spm_figure('GetWin','Figure 3'); clf    % display behavior
    spm_MDP_VB_game_tutoria_G_Copy(MDP); 

end