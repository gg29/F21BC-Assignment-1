%{
Authors: George Goniotakis - Murat Gunana
Course: F21BC Biologically Inspired Computation
Title: Part II - GA Version 3 COCO
Description: Includes HillClimbing, Mutation, Multiple Dimensions and
Crossover
%}

function xbest = GA_OPTIMIZER(FUN, DIM, ftarget, maxfunevals)

  popsize = 1000; %Population Size 
  maxfunevals = 500; %Maximum generations
  dimensions = DIM; %Dimensions
  crossoverRate = [0.2,0.8]; %[Rate of No Crossover, Rate of Crossover]
  mutationRate = [0.999,0.001]; %[Rate of No Mutation, Rate of Mutation]
  event = [0,1];
  fbest = inf; %Best Fitness so far
  mutationCounter = 0; %Number of Mutations
  crossoverCounter =0; %Number of Crossovers
  maxDigits = 5;
  
  individuals = randi([0 31], dimensions ,popsize);
  [rows, columns] = size(individuals);
    
  for currentGeneration = 1:ceil(maxfunevals/popsize)
             
    individuals = hillClimbing(individuals); %Do Hill-Climbing
    crossoverResult = randsample(event,1,true,crossoverRate); %Choose if mutation is going to happen
        
        %If a crossover is going to happen
        if crossoverResult == 1
           for r =1:rows:2
               crossoverBit = randi([2 maxDigits-1]); %Choose a random bit
              for c = 1:columns
                 element1 = de2bi(individuals(r,c), maxDigits, 'left-msb'); %Convert the vector to binary
                 element2 = de2bi(individuals(r+1,c), maxDigits, 'left-msb'); %Convert the vector to binary
                 crossPart1 = element1(crossoverBit:maxDigits);
                 crossPart2 = element2(crossoverBit:maxDigits);
                 element1(crossoverBit:maxDigits) = crossPart2;
                 element2(crossoverBit:maxDigits) = crossPart1;
                 individuals(r,c) = bi2de(element1, 'left-msb'); %Replace the new element to the table
                 individuals(r+1,c) = bi2de(element2, 'left-msb'); %Replace the new element to the table
              end
           end
           crossoverCounter = crossoverCounter +1; %Add 1 to the crossover counter
        end  
      
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

    