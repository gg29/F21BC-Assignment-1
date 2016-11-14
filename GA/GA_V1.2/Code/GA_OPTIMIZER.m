%{
Authors: George Goniotakis - Murat Gunana
Course: F21BC Biologically Inspired Computation
Title: Part II - GA Version 1.2 COCO
Description: Includes HillClimbing and Mutation
%}

function xbest = GA_OPTIMIZER(FUN, DIM, ftarget, maxfunevals)

  popsize = 1000; %Population size
  maxfunevals = 500; %Max Generations
  dimensions = DIM; %Dimensions
  mutationRate = [0.999,0.001]; %[Rate of No Mutation, Rate of Mutation]
  event = [0,1];
  fbest = inf; %Best Fitness so far
  mutationCounter = 0; %Number of Mutations
  maxDigits = 5; %Max Number of digits per individual
  
  individuals = randi([0 31], dimensions ,popsize);
  [rows, columns] = size(individuals);
    
  for currentGeneration = 1:ceil(maxfunevals/popsize)
             
    individuals = hillClimbing(individuals); %Do Hill-Climbing

    mutationResult = randsample(event,1,true,mutationRate); %Choose if mutation is going to happen
    
        %If a mutation is going to happen
        if mutationResult == 1
           for r =1:rows
              for c = 1:columns
                 mutationBit = randi(maxDigits); %Choose a random bit
                 element = de2bi(individuals(r,c), maxDigits, 'left-msb'); %Convert the vector to binary
                 element(:,mutationBit) = not(element(:,mutationBit)); %Flip this bit from 0 to 1 and vice versa
                 individuals(r,c) = bi2de(element, 'left-msb'); %Replace the new element to the table
              end
           end
           mutationCounter = mutationCounter +1; %Add 1 to the mutation counter
        end
      
    %xpop = 10 * rand(DIM, popsize) - 5;      % new solutions
    xpop = individuals;
    [fvalues, idx] = sort(feval(FUN, xpop)); % evaluate
  
    if fbest > fvalues(1)                    % keep best
      fbest = fvalues(1);
      xbest = xpop(:,idx(1));
    end
    if feval(FUN, 'fbest') < ftarget         % COCO-task achieved
      break;                                 % (works also for noisy functions)
    end
  end 
  
function [newMatrix] = hillClimbing(individuals)
    
    [rows, columns] = size(individuals);
    newMatrix = zeros(rows,columns);
    hillClimbingRate = 1;
    
    %Reduce each element by hillClimbingRate each time
    for r=1:rows
        for c=1:columns
            if individuals(r,c) ~= 0
                newMatrix(r,c) = individuals(r,c) - hillClimbingRate;
            end
        end
    end

    