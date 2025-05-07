function [X_rho_best,X_thick_best,OOBO_curve]=OOBOinvMT_func(InvParam,data)
% OOBOinvMT_func - Inversion of 1D Magnetotelluric (MT) data using the
% One-to-One-Based Optimizer (OOBO) algorithm.
%
% Inputs:
%   InvParam - Structure containing inversion parameters:
%              nPop   : Number of particles (population size)
%              niter  : Number of iterations
%              nlayer : Number of layers
%              rmin   : Minimum resistivity
%              rmax   : Maximum resistivity
%              tmin   : Minimum thickness
%              tmax   : Maximum thickness
%   data     - Observed MT data matrix [frequency, apparent resistivity, phase]
%
% Outputs:
%   X_rho_best   - Best estimated resistivity model (in log10 scale)
%   X_thick_best - Best estimated thickness model
%   OOBO_curve   - Convergence curve of best fitness value at each iteration

%% Preprocessing the observed data
freq =data(:,1); freq = freq(~isnan(freq)); %frequency;
app_data= data(:,2); app_data = app_data(~isnan(app_data)); % observed apparent resistivity
phase_data= data(:,3); phase_data = phase_data(~isnan(phase_data));%observed apparent phase
%% Initialization
% Initialize resistivity and thickness for all particles
X_rho=zeros(InvParam.nPop, InvParam.nlayer);
X_thick=zeros(InvParam.nPop, InvParam.nlayer-1);
for ipop = 1 : InvParam.nPop
    X_rho(ipop , :) = log10(InvParam.rmin) + rand*(log10(InvParam.rmax)- log10(InvParam.rmin));
    X_thick(ipop, :) = InvParam.tmin + rand*(InvParam.tmax- InvParam.tmin);
end
% Evaluate initial fitness for all particles
for ipop =1:InvParam.nPop
    L_rho=X_rho(ipop,:);
    L_thick=X_thick(ipop,:);
    [apparentResistivity, phase]=MT1D(10.^(L_rho),L_thick,freq);
    rhoapp_cal(ipop,:)=apparentResistivity;
    phase_mod(ipop,:)=phase;
    fitness=misfit(app_data,phase_data,rhoapp_cal(ipop,:),phase_mod(ipop,:));
    fit(ipop)=fitness;
end

%% Main OOBO Loop
for iter=1:InvParam.niter
    % Update the global best solution
    [best , location]=min(fit);
    if iter==1
        X_rho_best=X_rho(location,:);                   % Optimal location
        X_thick_best=X_thick(location,:);
        fbest=best;                    % The optimization objective function
    elseif best<fbest
        fbest=best;
        X_rho_best=X_rho(location,:);                                           % Optimal location
        X_thick_best=X_thick(location,:);
    end

    %% update agents position
    K=1:InvParam.nPop;
    for ipop=1:InvParam.nPop
        % Select a random peer k ≠ ipop
        k=randperm(InvParam.nPop-ipop+1,1);
        k=K(k);

        if k==ipop
            k=k+1;
            if k>InvParam.nPop
                k=k-2;
            end
        end
        K(find(K==k))=[]; % Remove selected index
        X_rho_k=X_rho(k,:);
        X_thick_k=X_thick(k,:);
        I=round(1+rand); % Binary value (1 or 2)

        % Position update equation based on relative fitness
        if fit(ipop)>fit(k)
            X_rho_new=X_rho(ipop,:)+ rand(1,InvParam.nlayer).*(X_rho_k-I.* X_rho(ipop,:));
            X_thick_new=X_thick(ipop,:)+ rand(1,InvParam.nlayer-1).*(X_thick_k-I.* X_thick(ipop,:));
        else
            X_rho_new=X_rho(ipop,:)+ rand(1,InvParam.nlayer).*(X_rho(ipop,:)-X_rho_k);
            X_thick_new=X_thick(ipop,:)+ rand(1,InvParam.nlayer-1).*(X_thick(ipop,:)-X_thick_k);
        end
        % Boundary handling for resistivity (log scale)
        U_rho=X_rho_new>log10(InvParam.rmax);        L_rho=X_rho_new<log10(InvParam.rmin);
        X_rho_new=(X_rho_new.*(~(U_rho+L_rho)))+log10(InvParam.rmax).*U_rho+log10(InvParam.rmin).*L_rho;
        % Boundary handling for thickness
        U_thick=X_thick_new>InvParam.tmax;        L_thick=X_thick_new<InvParam.tmin;
        X_thick_new=(X_thick_new.*(~(U_thick+L_thick)))+InvParam.tmax.*U_thick+InvParam.tmin.*L_thick;

         % Evaluate fitness of new solution
        [apparentResistivity, phase]=MT1D(10.^(X_rho_new),X_thick_new,freq);
        rhoapp_cal=apparentResistivity;
        phase_mod=phase;

        fitness=misfit(app_data,phase_data,rhoapp_cal,phase_mod);
        ff=fitness;

         % If fitness improved or equal, accept the new solution
        if ff <= fit (ipop)
            X_rho(ipop,:) = X_rho_new;
            X_thick(ipop,:) = X_thick_new;
            fit (ipop)=ff;
            %             S = struct('rho',X_rho, 'thick',X_thick,'error',fit,'misfit_curve',best_so_far);
            %             swarm(iter)=S;
        end
    end

    best_so_far(iter)=fbest;
    average(iter) = mean (fit);

    %     disp(['Iterasi ke-',num2str(iter)])
    % disp(['Progres : ', num2str((iter/InvParam.niter)*100), '%'])

end
Best_score=fbest;
OOBO_curve=best_so_far;

% save(sheetname,'swarm')

end