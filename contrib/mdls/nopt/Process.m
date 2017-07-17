%
% function Process(SC,DT,spar,mu,eta,overlap,testname)
%
% Name: Process
%
% Category: Auxiliar function
%
% Description: Process the data saved after the test and plots the required
% graphs. The tests are selected based in the input parameters that reduce
% the set of options.
%
% Input:
% SC ........ Sparse Coding algorithm selected. The options are: 'Lasso',
% 'MOL' or '*' (all)
% DT ........ Dictionary Training algorithm selected. The options are:
% 'ND', 'CG', 'MOCOD', the combination of both of them (i.e, 'NDCG') or '*'
% (all)
% spar ...... Array of values for sparsity.
% mu ........ Array of values for mu.
% eta ....... Array of values for eta.
% overlap ... Array of values for overlap.
% testname .. If a testname field is present the generated figures will be
% saved with <testname> as a prefix.
%
% Output:
% none
%
% Author: I. Ramirez and F. Lecumberry <nacho,fefo at fing dot edu dot uy>
%
% Version: $Idn$
%
function Process(SC,DT,spar,mu,eta,overlap,testname)
  home
  close all
  if isequal(nargin,5)
    saveimages = false;
  else
    saveimages = true;
  end
  %
  % Generated data files to load
  %  
  [datafiles,nF] =  ReadDataFiles();
  %
  % Generate minicode for printing
  %
  minicode = cell(nF,1);
  data = cell(nF,1);
  for k = 1:nF
    data{k} = load(datafiles{k});
    mc = ParseString(datafiles{k},'-');
    data{k}.overlapping = eval(mc{5}(3:4));
    data{k}.SC = mc{2};
    data{k}.DT = mc{3};
    data{k}.sparsity = eval(mc{6}(3:end));
    data{k}.mu = eval(mc{7}(3:end));
    data{k}.eta = eval(mc{8}(4:end));
    minicode{k} = strcat(data{k}.SC(1),'|',data{k}.DT,'|',...
      num2str(data{k}.sparsity),'|',...
      num2str(data{k}.mu),'|',...
      num2str(data{k}.eta));%,'|',num2str(data{k}.overlapping));
    data{k}.minicode = minicode{k};
    data{k}.fullcode = datafiles{k}(10:end-9); % EXTREMELY CHANCHO
  end
  %
  % Generate the list of valid experiments to label the graphs
  %
  ind = 1;
  for k = 1:nF
    if SelectTest(data{k},SC,DT,spar,mu,eta,overlap)
      fprintf(1,'%s\n',data{k}.minicode)
      legendArray{ind} = data{k}.minicode;
      ind = ind + 1;
    end
  end
  %
  % Define configuration for the graph 
  %
  mColor = 'rgbmkc';
  marker = '.sod*x+ph><^';
  markerSize = 5;
  lineWidth = 1;
  %
  % PSNR graph
  %
  hfp = Figure('PSNR');
  for k = 1:nF
    if SelectTest(data{k},SC,DT,spar,mu,eta,overlap)
      LC = length(mColor);
      col = mColor(mod(k-1,LC)+1);
      k2 = floor((k-1)/LC)+1;      
      mark = marker(mod(k2-1,length(marker))+1);
      plot(data{k}.psnrCurve,[col mark '-'],'MarkerSize',markerSize,'LineWidth',lineWidth)
      hold on
    end
  end
  legend(legendArray,'Location','SouthEast')
  ylabel('PSNR (dB)')
  xlabel('Iterations')
  title('PSNR')
  hold off
  grid on
  %
  % Coherence graph
  %
  hfcoh = Figure('Coherence');
  for k = 1:nF
    subplot(2,1,1)
    LC = length(mColor);
    col = mColor(mod(k-1,LC)+1);
    k2 = floor((k-1)/LC)+1;      
    mark = marker(mod(k2-1,length(marker))+1);
    if SelectTest(data{k},SC,DT,spar,mu,eta,overlap)
      plot(data{k}.cumCohCurve,[col mark '-'],'MarkerSize',markerSize,'LineWidth',lineWidth)
      hold on
    end
    subplot(2,1,2)
    if SelectTest(data{k},SC,DT,spar,mu,eta,overlap)
      plot(data{k}.maxCohCurve,[col mark '-'],'MarkerSize',markerSize,'LineWidth',lineWidth)
      hold on
    end
  end
  subplot(2,1,1)
  lhcoh = legend(legendArray,'Location','EastOutside');
  set(lhcoh, 'Position', [0.1 0.45 0.1 0.1])
  xlabel('Iterations')
  title('Cumulative Coherence')
  grid on
  hold off
  subplot(2,1,2)
  xlabel('Iterations')
  title('Maximal Coherence')
  grid on
  hold off
  %
  % Costs graphs
  %
  hfcst = Figure('Costs terms');
  for k = 1:nF
      LC = length(mColor);
      col = mColor(mod(k-1,LC)+1);
      k2 = floor((k-1)/LC)+1;      
      mark = marker(mod(k2-1,length(marker))+1);

    subplot(3,2,1)
    if SelectTest(data{k},SC,DT,spar,mu,eta,overlap)
      plot(data{k}.fittingCurve,[col mark '-'],'MarkerSize',markerSize,'LineWidth',lineWidth)
      hold on
    end
    subplot(3,2,2)
    if SelectTest(data{k},SC,DT,spar,mu,eta,overlap)
      plot(data{k}.regCurve,[col mark '-'],'MarkerSize',markerSize,'LineWidth',lineWidth)
      hold on
    end
    subplot(3,2,3)
    if SelectTest(data{k},SC,DT,spar,mu,eta,overlap)
      plot(data{k}.cohCurve,[col mark '-'],'MarkerSize',markerSize,'LineWidth',lineWidth)
      hold on
    end
    subplot(3,2,4)
    if SelectTest(data{k},SC,DT,spar,mu,eta,overlap)
      plot(data{k}.normCurve,[col mark '-'],'MarkerSize',markerSize,'LineWidth',lineWidth)
      hold on
    end
    subplot(3,2,[5 6])
    if SelectTest(data{k},SC,DT,spar,mu,eta,overlap)
      plot(data{k}.costCurve,[col mark '-'],'MarkerSize',markerSize,'LineWidth',lineWidth)
      hold on
    end
  end
  subplot(3,2,1)
  xlabel('Iterations')
  title('Fitting term')
  grid on
  hold off
  subplot(3,2,2)
  % legend(legendArray,'Location','EastOutside')
  xlabel('Iterations')
  title('Regularization term')
  grid on
  hold off
  subplot(3,2,3)
  % legend(legendArray,'Location','EastOutside')
  xlabel('Iterations')
  title('Coherence term')
  grid on
  hold off
  subplot(3,2,4)
  % legend(legendArray,'Location','EastOutside')
  xlabel('Iterations')
  title('Unnormalization term')
  grid on
  hold off
  subplot(3,2,[5 6])
  % legend(legendArray,'Location','EastOutside')
  legend(legendArray,'Location','EastOutside');
  %   set(lhcst, 'Position', [0.01 0.5 0.1 0.1])
  xlabel('Iterations')
  title('Cost term')
  grid on
  hold off
  
  hfcst = Figure('Dictionary evolution');
  for k = 1:nF
      LC = length(mColor);
      col = mColor(mod(k-1,LC)+1);
      k2 = floor((k-1)/LC)+1;      
      mark = marker(mod(k2-1,length(marker))+1);
    if SelectTest(data{k},SC,DT,spar,mu,eta,overlap)
      codename = data{k}.fullcode;
      iters = 1;
      while true
        fname = sprintf('../output/Dictionaries/%s-Dictionary-i%02d.mat',codename,iters);
        if ~exist(fname,'file')
          iters = iters - 1;
          break;
        end
        iters = iters + 1;
      end
      if iters > 0                
        fname = sprintf('../output/Dictionaries/%s-Dictionary-i%02d.mat',codename,iters);
        Daster = load(fname);      
        err = zeros(1,iters);
        for i = 1:iters
          fname = sprintf('../output/Dictionaries/%s-Dictionary-i%02d.mat',codename,i);
          Di = load(fname);
          err(i) = norm(Daster.D-Di.D,'fro');
        end
        plot(err,[col mark '-'],'MarkerSize',markerSize,'LineWidth',lineWidth) 
      else
        fprintf(1,'Could not find dictionary sequence for %s\n',codename);
        hold on
        plot([0],[col mark '-'],'MarkerSize',markerSize,'LineWidth',lineWidth); 
      end
    end  
  end
  xlabel('Iteration');
  ylabel('||D-Dfinal||');
  legend(legendArray,'Location','EastOutside');

  if saveimages
    saveas(hfp,sprintf('../output/%s-psnrs.png',testname))
    saveas(hfp,sprintf('../output/%s-psnrs.fig',testname))
    saveas(hfcoh,sprintf('../output/%s-coherences.png',testname))
    saveas(hfcoh,sprintf('../output/%s-coherences.fig',testname))
    saveas(hfcst,sprintf('../output/%s-costs.png',testname))
    saveas(hfcst,sprintf('../output/%s-costs.fig',testname))
    %   saveas(hfn,sprintf('../output/%s-norms.png',testname))
    %   saveas(hfn,sprintf('../output/%s-norms.fig',testname))
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [datafiles,nF] =  ReadDataFiles()
  files = dir('../output/*-data.mat');
  nF = length(files);
  datafiles = cell(nF,1);
  for k = 1:nF
    datafiles{k} = strcat('../output/',files(k).name);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out =  SelectTest(data,SC,DT,spar,mu,eta,overlap)
  %
  % Text fields
  %
  if ~isempty(strfind(SC,data.SC)) || isequal(SC,'*')
    t1 = true; else t1 = false;
  end
  if ~isempty(strfind(DT,data.DT)) || isequal(DT,'*')
    t2 = true; else t2 = false;
  end
  %
  % Numeric fields
  %
  if isequal(spar,'*') || any(spar == data.sparsity)
    t3 = true; else t3 = false;
  end
  if isequal(mu,'*') || any(mu == data.mu)
    t4 = true; else t4 = false;
  end
  if isequal(eta,'*') || any(eta == data.eta)
    t5 = true; else t5 = false;
  end
  if isequal(overlap,'*') || any(overlap == data.overlapping)
    t6 = true; else t6 = false;
  end
  %
  % Boolean output
  %
  out = t1 & t2 & t3 & t4 & t5 & t6;
end
