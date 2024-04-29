%% Trial to model an active inference agent in a maze with two hints

clear all
close all      % These commands clear the workspace and close any figures

rng('shuffle') % This sets the random number generator to produce a different 
               % random sequence each time, which leads to variability in 
               % repeated simulation results (you can alse set to 'default'
               % to produce the same random sequence each time)

%% SETTING UP

% Number of time points or 'epochs' within a trial: T
% =========================================================================

% Here, we specify 4 time points (T), in which the agent 1) starts in a 'Start'
% state, 2) first moves to either the first hint or goes direclty for
% Left/Right, 3) moves from the 1st hint to the 2nd hint or from the 1st hint
% to Left/Right, or from the Left/Right to the start state, 4) either moves
% from the 2nd hint to Left/Right or from Left/Right back to start. 

T = 4;


rs1 = 4; % Risk-seeking parameter (set to the variable rs below)


%==========================================================================
% Set up the state factors D and d
%==========================================================================
% 'Context' state factor

D{1} = [1, 0]'; % {'left better', 'right better'}

% 'Behaviour' state factor the agent always begins in the start state

D{2} = [1 0 0 0 0]'; %{'Start', 'get 1st hint', 'get  2nd hint', 'Left', 'Right'}

% These are transposed so that they become column vectors and this is why
% later you find left-better and right-better in columns and the same goes
% for the actions (check B)

% Prior beliefs over intial states (d). This is optional and generates
% learning priors over states if specified. 

% For context beliefs the agents starts by believing that both states are
% equally probable but with low confidence. 

d{1} = [.25, .25]'; % {'left better','right better'}

% For behaviour beliefs the agent is certain that it will begin in the
% start state 

d{2} = [1 0 0 0 0 ]'; %{'Start', 'get 1st hint', 'get 2nd Hint', 'Left', 'Right'}

%==========================================================================
% State Outomes Mapping and Beliefs A and a
%==========================================================================

%Ns = [length(D{1}) length(D{2})]

% I start by specifying that both context D{1} generate the 'No Hint'
% observation across all behaviour states (1, 2, 3, 4, 5) 

A{1}(:, :, 1) =[1 1;   % No hint                
                0 0;   % get 2nd hint
                0 0;   % Left
                0 0];  % Right
               
A{1}(:, :, 2) =[1 1;   % No hint                
                0 0;   % get 2nd hint
                0 0;   % Left
                0 0];  % Right

A{1}(:, :, 3) =[1 1;   % No hint                
                0 0;   % get 2nd hint
                0 0;   % Left
                0 0];  % Right

A{1}(:, :, 4) =[1 1;   % No hint                
                0 0;   % get 2nd hint
                0 0;   % Left
                0 0];  % Right

A{1}(:, :, 5) =[1 1;   % No hint                
                0 0;   % get 2nd hint
                0 0;   % Left
                0 0];  % Right




%--------------------------------------------------------------------------
% Then I specify that the 'Get Hint' behavior (2) state generates a hint that
% 50% it can take the 2nd hint in both context, or it can generate a go
% left or go right  depending on the context state. In the first hint case, 
% the hints are accurate with a probability of pHA. 


pHA = 1; 


A{1}(:, :, 2) = [0          0;   % No hint                
                 pHA/2  pHA/2;   % get 2nd hint
                 pHA   1-pHA;   % Left
                 1-pHA  pHA];  % Right

% Then I specify the get 2nd hint behaviour, which goes with the chance prob 
% but has lower probability of telling the truth.

pHB = 0.7;

A{1}(:, :, 3) =  [0         0;  % No hint
                  0         0;  % get 2nd Hint 
                  pHB 1-pHB;  % Left 
                  1-pHB pHB]; % Right

%--------------------------------------------------------------------------
% Next I specify the win and losses. First the three behaviour 'start' and 'get 1st hint'
% and get '2nd hint' do not map to any osservation of win and losses 

A{2}(:, :, 1) = [1 1;   % Null
                0 0;   % Punishment
                0 0];  % Reward

A{2}(:, :, 2) = [1 1;   % Null
                0 0;   % Punishment
                0 0];  % Reward

A{2}(:, :, 3) = [1 1;   % Null
                0 0;   % Punishment
                0 0];  % Reward

% Now I specify how the behaviour 'choose left' and 'choose right'
% generates the loss or win with a probability of pWin.

pWin = 1; 

% Here the agent gets punishment or reward with 100% prob depending on the
% context


A{2}(:, :, 4) = [0           0;    % Null
                 1-pWin pWin;    % Punishment
                 pWin 1-pWin];   % Reward

% Here the agent 'choose right' (4), hence the probabilities of reward are
% inverted.

A{2}(:, :, 5) = [0          0;     % Null
                pWin 1-pWin;     % Punishment
                1-pWin pWin];    % Reward

%--------------------------------------------------------------------------
% Finally the agent knows that behaviours were carried on as planned, that
% is, it observes its own behaviours. It basically knows that when it is in
% the start state it is in the start state and so on. 

A{3}(:, :, 1) = [1 1;  % Start
                 0 0;  % Get 1st hint
                 0 0;  % Get 2nd Hint
                 0 0;  % Choose left
                 0 0]; % Choose right   

A{3}(:, :, 2) = [0 0;  % Start
                 1 1;  % Get 1st hint
                 0 0;  % Get 2nd Hint
                 0 0;  % Choose left
                 0 0]; % Choose right 

A{3}(:, :, 3) = [0 0;  % Start
                 0 0;  % Get 1st hint
                 1 1;  % Get 2nd Hint
                 0 0;  % Choose left
                 0 0]; % Choose right 

A{3}(:, :, 4) = [0 0;  % Start
                 0 0;  % Get 1st hint
                 0 0;  % Get 2nd Hint
                 1 1;  % Choose left
                 0 0]; % Choose right 

A{3}(:, :, 5) = [0 0;  % Start
                 0 0;  % Get 1st hint
                 0 0;  % Get 2nd Hint
                 0 0;  % Choose left
                 1 1]; % Choose right 

               
%--------------------------------------------------------------------------

%==========================================================================
% Here I can insert a{} to specify prior beliefs about state-outcome mapping
% in the generative model. (a) This is optional and if specified will
% simulate learning state-outcome mapping, for exmample, learning the reward
% probabilities or learning the hint accuracy. 
%==========================================================================

%==========================================================================
% Controlled transitions B{} and transitions beliefs b{}
%==========================================================================
% Here I specify the probabilitisc transitions between hidden states under
% each action. Because there are two state factors D{1} and D{2}, I need two
% sets of matrices, one for each state factor.

% Columns are state at time t, rows are states at time t+1

% The agent cannot control the context state factor, so there is only 1 'action'
% indicating that the context remains stable within a trial, in other words
% states at time t (columns) remain the same at state t + 1 (rows). There
% is only one 'action' possible, that is transitioning from one state to
% the other, hence the third dimensions does not change, there is no 
% B{1}(:,:,2) etc. 

B{1}(:,:,1) = [1, 0;  % Left-better
               0, 1]; % Right-better
%--------------------------------------------------------------------------

% The agent can instead control the behaviour state factor, which has 5
% possible actions at each time step. Here therefore I need 5 matrices.
% Notice how it is described as moving TO the state instead that FROM.
% NOTE: EACH COLUMN MUST CONTAIN A 1, IT'S MANDATORY

% Move to the start state from any other state 

B{2}(:, :, 1) = [1 1 1 1 1;   % Start state
                 0 0 0 0 0;   % get 1st Hint
                 0 0 0 0 0;   % get 2nd Hint 
                 0 0 0 0 0;   % Choose Left
                 0 0 0 0 0];  % Choose Right

% Move to the 1st Hint state from any other state but from the 2nd hint

B{2}(:, :, 2) = [0 0 0 0 0;   % Start state
                 1 1 1 1 1;   % get 1st Hint
                 0 0 0 0 0;   % get 2nd Hint 
                 0 0 0 0 0;   % Choose Left
                 0 0 0 0 0];  % Choose Right

% Move to the take 2nd hint from any other state but from start state and
% from choose left and right, bc the location of the second hint is
% unveiled by the first hint

B{2}(:, :, 3) = [0 0 0 0 0;   % Start state
                 0 0 0 0 0;   % get 1st Hint
                 1 1 1 1 1;   % get 2nd Hint 
                 0 0 0 0 0;   % Choose Left
                 0 0 0 0 0];  % Choose Right

% Move to the Choose Left state from any other state 

B{2}(:, :, 4) = [0 0 0 0 0;   % Start state
                 0 0 0 0 0;   % get 1st Hint
                 0 0 0 0 0;   % get 2nd Hint 
                 1 1 1 1 1;   % Choose Left
                 0 0 0 0 0];  % Choose Right

% Move to the Choose Right State from any other state

B{2}(:, :, 5) = [0 0 0 0 0;   % Start state
                 0 0 0 0 0;   % get 1st Hint
                 0 0 0 0 0;   % get 2nd Hint 
                 0 0 0 0 0;   % Choose Left
                 1 1 1 1 1];  % Choose Right

%==========================================================================
% Again I could here specify b{} that is prior beliefs about state
% transitions in the generative model (b). This is a set of matrices with
% the same structure as B. This is optional and again will simulate
% learning state transitions if specified. 
%==========================================================================

%==========================================================================
% Preferred outcomes C and c
%==========================================================================
% Here I specify preferences over each set of outcomes (C), with one matrix
% per outcome modality (numbers of A's). Here columns indicate time points, and rows
% indicate observations (same order as in the corresponding A matrices).
% Therefore columns ALWAYS = T

% I can start by setting a 0 preference for all outcomes: 

% No = [size(A{1},1) size(A{2},1) size(A{3},1)], this would take the number
% of outomes per modality, if printed returns (3 , 3 , 5) which basically
% means that the first A{1} and A{2} have 3 rows (number of outcomes) and
% A{3} has 5 rows (number of outcomes). The matrices basically are (number
% of outcomes/Rows, time points/Columns)
% hence C{1} = (3, 3), C{2} = (3, 3), C{3} = (5, 3)
% 
% One could then put 
%C{1}     = zeros(No(1),T); % Hints
%C{2}      = zeros(No(2),T); % Wins/Losses
%C{3}      = zeros(No(3),T); % Observed Behaviors
% I will set them manually here

C{1} = [0 0 0 0;  % No Hint
        0 0 0 0;  % get 2nd hint
        0 0 0 0;  % Left Hint
        0 0 0 0]; % Right Hint


C{2} = [0 0 0 0;  % Null
        0 0 0 0;  % Punishment
        0 0 0 0]; % Reward

C{3} = [0 0 0 0;  % Start
        0 0 0 0;  % Get 1st hint
        0 0 0 0;  % Get 2nd Hint
        0 0 0 0;  % Choose left
        0 0 0 0]; % Choose right 

%--------------------------------------------------------------------------

% Now after setting them to zeros one can decide whether to change them, in
% the sense whether to insert preferences over outcomes. The question would be
% "Does the agent have a preference over No Hint, Left Hint, Right Hint?" 
% "Does the agent have a preference over Null, Punishment, and Reward?"
% "Does the agent have a preference over "Start, Get 1st hint, get 2nd hint C. Left, C.
% Right?"
% Here I will only specify that the agent prefer to avoid punishiment and
% seeks reward, however one could for example, think that if a person is
% left-handed there might be some preference for choosing left vs choosing
% right hence could specify also C{1} and C{3}.

% I can specify a punishmentaversion magnitude (pa) and a risk seeking magnitude
% or 'reward seeking' % (rs), at time poits 2 and 3. What would happen if
% having this set also at time point 1? 

rs = rs1; % I set this at the top of the script. Higher values promote risk seeking

pa = 1; 

% The first row goes all at 0 as it econde the preference over observing a
% null, which is not what is interesting. The second row encodes punishment
% hence will be -pa to signify a preference against observing a loss at
% times 2 and 4. The third row, is going to be rs, signifying a preference
% versus observing a reward at times 2 and 3

C{2} = [0   0   0    0;  % Null
        0 -pa -pa  -pa;  % Punishment
        0  rs  rs   rs]; % Reward

%==========================================================================
% Again I could specify a c{}, that is to simulate preference learning,
% specifically through a Dirichlet distribution over preferences (c).
% This works by increasing the preference magnitued for an outcome each
% time that outcome is observed. The assumption is that preferences
% naturally increase forentering situations that are more familiar. 
%==========================================================================

%==========================================================================
% Shallow vs Deep Policies
%==========================================================================

% Policies are encoded in vectors, the number of vectors is determined by
% the number of state factors D. In this task I inserted four time points,
% T = 4 which means that a policy will consist of two actions. 

% Shallow policies are policies where the model looks only one step ahead
% and are econded in the vector U, where we have one column for each allowable
% action. I am not going to use the vector U in
% this example, however, it would look like this: 

% U(:, :, 1) = [1 1 1 1]; % Context is not controllable 
% U(:, :, 2) = [1 2 3 4]; % All actions are allowed

% I could also specify it like:
% Np = 4; % Number of policies
% Nf = 2; % Number of state factors
% U = ones(1, Np, Nf);


% Deep Policies are policies in which the model looks ahead until the end
% of the trial. Hence I need to specify one column for each allowable
% policy (with each entry inidicating an action number) with one row per
% time point and a third dimension specifying each state factor. 
% In other words, rows correspond to time points and should be length T-1
% (here, 3 transitions, from time point 1 % to time point 2, 
% fomr time point 2 to time point 3, and from 3 to 4). ALWAYS ONE ROW LESS
% BC IS THE ACTIONS THAT YOU DON'T DO AT THE FRIST STEP


V(:, :, 1) = [1 1 1 1 1 1 1;
              1 1 1 1 1 1 1; 
              1 1 1 1 1 1 1]; % Context is not controllable

V(:, :, 2) = [1 2 2 2 2 1 1;
              1 3 3 4 5 4 5;
              1 4 5 1 1 1 1];

% For V(:,:,2), columns left to right indicate policies allowing: 
% 1. staying in the start state 
% 2. taking the 1st hint then the 2nd hint then Left 
% 3. taking the 1st hint then the 2nd hint then Right
% 4. taking the 1st hint then Left then back 
% 5. Takin ghte 1st gint the right then back 
% 6. Choosing the left and then back 
% 7. Choosing the right and then back 




% I could also specify it as: 
% Np = 5; % Number of policies as higlighted above
% Nf = 2; % Number of state factors
% V = ones(T-1,Np,Nf)

%==========================================================================
% Fixed priors overe policies/Habits: E and e. 
%==========================================================================

%--------------------------------------------------------------------------
% Optional: a columns vector with one entry per policy, indicating the 
% prior probability of choosing that policy (i.e., independent of other 
% beliefs). 
%--------------------------------------------------------------------------

% We will not equip our agent with habits in our example simulations, 
% but this could be specified as a follows if one wanted to include a
% strong habit to choose the 4th policy:

% E = [.1 .1 .1 .6 .1]';

% To incorporate habit learning, where policies become more likely after 
% each time they are chosen, one can also specify concentration parameters
% by specifying e. For example:

% e = [1 1 1 1 1]';

%% Define MDP Structure 
%==========================================================================

mdp.t = T; % Number of time steps

mdp.V = V; % Allowable deep policies

mdp.A = A; % state-outcome mapping

mdp.B = B; % Transition probabilities

mdp.C = C; % preferred states 

mdp.D = D; % prior over initial states

mdp.d = d; % enable learning priors over intial states

%--------------------------------------------------------------------------
% Labelling and plotting

% Adding some lables for subsequent plotting 

label.factor{1} = 'Context'; label.name{1} = {'Left Better', 'Right Better'};
label.factor{2} = 'Choiche States'; label.name{2} = {'Start', 'Get 1st Hint', 'Get 2nd Hint', 'Choose Left', 'Choose Right'};
label.modality{1} = 'Hint'; label.outcome{1} = {'No hint', 'Get 2nd Hint', 'Choose Left', 'Choose Right'};
label.modality{2} = 'Reward/Punishment'; label.outcome{2} = {'Null', 'Punishment', 'Reward'};
label.modality{3} = 'Observed Action'; label.outcome{3} = {'Start', 'Get 1st Hint', 'Get 2nd Hint', 'Choose Left', 'Choose Right'};
label.action{2} = {'Start', '1st Hint', '2nd Hint', 'Left', 'Right'};
mdp.label = label;


% Check whether all matrix dimensions are correct 

mdp = spm_MDP_check(mdp);

% Run active inference estimantion 

MDP = spm_MDP_VB_X(mdp);

% We can then use standard plotting routines to visualize simulated neural 
% responses

spm_figure('GetWin','Figure 1'); clf    % display behavior
spm_MDP_VB_LFP(MDP); 

%  and to show posterior beliefs and behavior:

spm_figure('GetWin','Figure 2'); clf    % display behavior
spm_MDP_VB_trial(MDP); 























