%{
Authors: George Goniotakis - Murat Gunana
Course: F21BC Biologically Inspired Computation
Title: Part II - GA Version 1
Description: Includes HillClimbing and Mutation
%}

function Fvalue = GA_V1()

    clc,clear,close all; %Clear command window
    maxDigits = 10; %Maximum Number of Binary Digits
    popsize = 100; %Population Size 
    generations = 100; %Maximum Number of Generations
    dimensions = 2; %Dimensions
    mutationRate = [0.999,0.001]; %[Rate of No Mutation, Rate of Mutation]
    event = [0,1];
    bestFitness = inf; %Best Fitness so far
    ftarget = 0; %Target Value
    success = 0; %If global value achieved
    mutationCounter = 0; %Number of Mutations
    currentGeneration = 0; %Generations Counter

    %20 Different iterations to find random numbers
    for i=1:20
        indBinary = randi([0 1], popsize*dimensions, maxDigits); %Create a binary array with 2 individuals
    end

    individuals= bi2de(indBinary, 'left-msb'); %The population in binary is converted to decimal
    individuals = reshape(individuals, popsize, dimensions); %Reform the vector
    [rows, columns] = size(individuals);
    
    %Calculate the fitness for each pair for the first time
    for r= 1: rows
        tempScore = sphereFunction(individuals(r,:));
        
        if tempScore < bestFitness
            bestFitness = tempScore;
        end
    end

    %Start the hillClimbing algorithm until global value is achieved
    while currentGeneration < generations && success == 0 
        
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
        
        %Evaluate the new table
        for r= 1: rows
           tempScore = sphereFunction(individuals(r,:));

           if tempScore < bestFitness
              bestFitness = tempScore;
           end
        end
        
        %Check if global value has achieved
        if bestFitness == ftarget
           success = 1;
           break;
        end

        currentGeneration = currentGeneration +1; %Increase the generation counter
    end 
    
    %If the check is successful
    if success == 1
        disp(sprintf(['Best fitness found after %d generations and %d mutations'],...
                   currentGeneration,mutationCounter));
    %Else if total generations exceeded
    else
        disp(sprintf(['Total available generations exceeded. Best fitness so far %e. Mutations done %d']...
            ,bestFitness,mutationCounter));
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
end
