%{
Authors: Murat Gunana - George Goniotakis
Course: F21BC Biologically Inspired Computation
Title: Part I – MLP_GA Final Version COCO
Description: Multiple Dimensions, Crossover, Roulette Selection
%}

function xbest = MLP_GA_OPTIMIZER(FUN, DIM, ftarget, maxfunevals)
% MY_OPTIMIZER(FUN, DIM, ftarget, maxfunevals)
% samples new points uniformly randomly in [-5,5]^DIM
% and evaluates them on FUN until ftarget of maxfunevals
% is reached, or until 1e8 * DIM fevals are conducted. 
  
  maxfunevals = 200; 
  fbest = inf;
   
  % Builds ANN topology.
  inputs = 'Please enter number of inputs: ';
  input_no = input(inputs);
  layer = 'Please enter number of layers:  ';
  hidden_layer = input(layer);
  
  % Constructs possible number of genes for output based on inputs. 
  outputVector = 2^input_no;
  desired = sprintf('Please enter the desired output (between 0 to %d)',(outputVector -1));
  desired_output = input(desired);
    
  n = 'Please enter the population size:  ';
  popsize = input(n);
  
  output = 1;
  weightV = (input_no + 1) * (input_no*hidden_layer + output);
  
  % MLP-GA
  MAX_INPUT = (2.^weightV - 1);
  phen = randi([0 MAX_INPUT-1], DIM, popsize);
   
  % Calculate the best overall fitness of the population.
  totalFitness = 0;
  bestGen = 1;
  for iter = 1:ceil(maxfunevals/popsize)
    % This bit of code selects the individuals based on roulette rule
    % selection with their fitness values assigned. The rule randomly 
    % select each individual with probabilites assigned with fitness
    % function. It then pairs the fittest individuals to perform crossover.
    [rows, columns] = size(phen);
    %crossMatrix = randi([0 0],1, rows*columns);
    populationVector = reshape(phen, 1, rows*columns);
    
    % Reproduction based on roulette selection.
    crossMatrix = rouletteSelection(populationVector);
    reproductionVector = zeros(rows, columns);
    for rSize = 1 : columns
      if length(crossMatrix(1,:)) > 1
        element = randsample(crossMatrix, 1);
        crossMatrix = indSelection(element, crossMatrix);
      else
        element = crossMatrix;
      end
      reproductionVector(1,rSize) = element;
    end
 
    % Crossover operations.
    [rows, columns] = size(reproductionVector);
    for xSize = 0 : 2 : columns - 2
      event = [0,1];
      crossoverRate = [0.4,0.6];
      crossoverResult = randsample(event,1,true,crossoverRate);
      if crossoverResult == 1
        parent1 = de2bi(reproductionVector(1,(xSize+1)), weightV);
        parent2 = de2bi(reproductionVector(1,(xSize+2)), weightV);
        
        rNumber = randi([3 (weightV-2)]);
        crossPoint = rNumber(1,1);
        crossPointLenght = weightV - crossPoint;
        
        sizel = crossPointLenght+1;
        xover1 = randi([0 0],1, sizel);
        xover2 = randi([0 0],1, sizel);
        xover1(1:sizel) = parent1(crossPoint:end);
        xover2(1:sizel) = parent2(crossPoint:end);
 
        parent1(crossPoint:end) = [];
        parent2(crossPoint:end) = [];
        child1 = horzcat(parent1,xover2);
        child2 = horzcat(parent2,xover1);
 
        reproductionVector(1,xSize+1)= bi2de(child1, 'left-msb');
        reproductionVector(1,xSize+2)= bi2de(child2, 'left-msb');
      end
    end
 
    % Mutation operation happens here. Ideally the mutaion should happen in a 
    % fixed point rather than randomly selected. In this example it happens on
    % third bit but can be assigned to a random number which in that case is
    % the "mutate" variable.
    for row = 1 : rows
      for column = 1 : columns
        mutationRate = [0.999,0.001];
        event = [0,1];
        mutationResultResult = randsample(event,1,true,mutationRate);
        if mutationResultResult == 1
          element = de2bi(reproductionVector(row,column), weightV);
          mutate = randi(weightV);
          element(:,mutate) = not(element(:,mutate));
          reproductionVector(row, column) = bi2de(element, 'left-msb');
        end
      end
    end
    if sum(reproductionVector.^2) > totalFitness
      totalFitness = sum(reproductionVector.^2);
      bestGenPop = reproductionVector;
      bestGen = i;
    end
    % Gets the best individual chosen from best generation. At this point we
    % are certain we found the best weights for our neural network. Now we
    % can calculate overall error between the best individual output which 
    % we tried to maximize the function as 2^input_no. We can test the
    % desired output against the neural network which will return the error
    % rate.
    [rows,columns] = size(reproductionVector);
    for row = 1: rows
      for column = 1: columns
        [delta] = NN_config(reproductionVector(row,column), input_no, output, desired_output, hidden_layer, weightV);
        desiredMatrix = zeros(length(delta(:,1)),1);
        for deltaSize = 1: length(delta(:,1))
          if delta(deltaSize,1) < 0
            desiredMatrix(deltaSize,1) = 0;
          else
            desiredMatrix(deltaSize,1) = 1;
          end
        end
      end
    end
    
   % xpop = 10 * rand(DIM, popsize) - 5;      % new solutions
   xpop = reproductionVector;
    [fvalues, idx] = sort(feval(FUN, xpop)); % evaluate
    if fbest > fvalues(1)                    % keep best
      fbest = fvalues(1);
      xbest = xpop(:,idx(1));
    end
    if feval(FUN, 'fbest') < ftarget         % COCO-task achieved
      break;                                 % (works also for noisy functions)
    end
  end 
  
  % Set up neural network configurations.
function [delta] = NN_config(individual, input_no, output_neuron, desired_output, hidden_layer, weightV);
  % Builds input matrix.
  input = dec2bin(0:(2^input_no-1));
  % Initialize the bias
  bias = ones(input_no*hidden_layer+output_neuron,1)*(1);
  desired_output = de2bi(desired_output, 2^input_no);
  % Encodes the weights.
  chrom = de2bi(individual, weightV);
  chrom = reshape(chrom, (input_no*hidden_layer + output_neuron), (input_no +1));
 
  % Converts weight bits to integers.
  weights = randi([ 0 0], (input_no*hidden_layer + output_neuron), (input_no + 1));
  weights = -1*2.*weights;
  [rows, columns] = size(chrom);
  for i = 1: rows
    for j = 1: columns
      if chrom(i, j) == 1
        weights(i, j) = 1;
      else
        weights(i, j) = -1;
      end
    end
  end
  
  out = zeros(2^input_no,1);
  delta = zeros(2^input_no,1);
  inputLength = length(input(:,1));

  for j = 1:inputLength
    
    for k = 1: hidden_layer + 1
      layer_input = zeros(input_no, 1);
      if k == 1
        % If it's the first layer it should take the input directly,
        % otherwise it will use previous layer neuron uotputs.
        for i = 1: input_no
          layer_input(i, 1) = bin2dec(input(j,i));
        end
      end
      
      % Ensures that it's the output layer.
      if k == hidden_layer + 1 
        layerNeuronNumber = 1;
      else
        layerNeuronNumber = input_no;
      end
      
      % For each neuron in the hidden layer, it calculates the output.
      for l = 1: layerNeuronNumber
        layer_neuron_output = 0;
        for m = 1: input_no + 1
          if m == 1
            % Bias case.
            layer_neuron_output = layer_neuron_output + bias(l,m)*weights(m,1);
          else
            % For input or output of a neuron.
            layer_neuron_output = layer_neuron_output + layer_input(m-1,1)*weights(l,m);
          end
        end
        layer_input(l,1) = 1/(1+exp(-layer_neuron_output));
        out(j,1) = 1/(1+exp(-layer_neuron_output));
        delta(j,1) = out(j,1)*(1-out(j,1))*(desired_output(1,j)-out(j,1));
      end
    end
  end
  
% Eliminates the element from the selection matrix and returns the updated
% matrix.
function matrix = indSelection(element, matrix)
  [rows, columns] = size(matrix);
  for row = 1 : rows
    for column = 1 : columns
      if matrix(row,column) == element
        matrix(:,column) = [];
        break;
      end
    end
  end
 
% Calculates the fitness value of each individual and makes a random 
% selection based on the wieghts which consist of their fitness value. 
function crossMatrix = rouletteSelection(population)
  [rows, columns] = size(population);
  crossMatrix = zeros(rows, columns);
  sumFitness = sum(population);
  weight = population./sumFitness;
 
  for row = 1: rows
    for column = 1: columns
      ind = randsample(population,1,true,weight);
      crossMatrix(row,column) = ind;
    end
  end
 
  


