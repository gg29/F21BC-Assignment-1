% runs an entire experiment for benchmarking MY_OPTIMIZER
% on the noise-free testbed. fgeneric.m and benchmarks.m
% must be in the path of Matlab/Octave
% CAPITALIZATION indicates code adaptations to be made

clc;
addpath('C:\Users\Dell\Documents\MATLAB\MLP\Code');  % should point to fgeneric.m etc.
datapath = 'C:\Users\Dell\Documents\MATLAB\MLP\Code\experiments\exp1';  % different folder for each experiment
% opt.inputFormat = 'row';
opt.algName = 'Sphere Function';
opt.comments = '';
maxfunevals = '10 * dim'; % 10*dim is a short test-experiment taking a few minutes 
                          % INCREMENT maxfunevals successively to larger value(s)
minfunevals = 'dim + 2';  % PUT MINIMAL SENSIBLE NUMBER OF EVALUATIONS for a restart
maxrestarts = 1e4;        % SET to zero for an entirely deterministic algorithm

dimensions = [2,3,5,10,20,40];  % small dimensions first, for CPU reasons
functions = 1; %benchmarks('FunctionIndices');  % or benchmarksnoisy(...)
instances = [1:5, 41:50];  % 15 function instances
averageError = 0;

more off;  % in octave pagination is on by default

t0 = clock;
rand('state', sum(100 * t0));

for dim = dimensions
  for ifun = functions
    for iinstance = instances
      fgeneric('initialize', ifun, iinstance, datapath, opt); 

      % independent restarts until maxfunevals or ftarget is reached
      for restarts = 0:maxrestarts
        if restarts > 0  % write additional restarted info
          fgeneric('restart', 'independent restart')
        end
        MLP_GA_OPTIMIZER('fgeneric', dim, fgeneric('ftarget'), ...
                     eval(maxfunevals) - fgeneric('evaluations'));
        if fgeneric('fbest') < fgeneric('ftarget') || ...
           fgeneric('evaluations') + eval(minfunevals) > eval(maxfunevals)
          break;
        end  
      end

       disp(sprintf(['  f%d in %d-D, instance %d: FEs=%d with %d restarts,' ...
                     ' fbest-ftarget=%.4e, elapsed time [h]: %.2f'], ...
                    ifun, dim, iinstance, ...
                    fgeneric('evaluations'), ...
                    restarts, ...
                    fgeneric('fbest') - fgeneric('ftarget'), ...
                    etime(clock, t0)/60/60));
        averageError = averageError + fgeneric('fbest') - fgeneric('ftarget');
      fgeneric('finalize');
    end
    
    averageError = averageError/length(instances);
    disp(sprintf([' Average Error for f%d in %d-D is %d'],...
                   ifun, dim, averageError));
     disp(['      date and time: ' num2str(clock, ' %.0f')]);
  end
   disp(sprintf('---- dimension %d-D done ----', dim));
end

