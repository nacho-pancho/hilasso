function CreateNoisyImages(listOfImages,listOfSigmas)
  if ~exist('listOfSigmas','var')
    listOfSigmas = [ 10 ];
  end
  if ~exist('listOfImages','var')
    listOfImages = { ...
      '../images/barbara.png',...
      '../images/boat.png',...
      '../images/mandrill.png',...
      '../images/man.png',...
      '../images/lena.png',...
      '../images/peppers.png',...
      '../images/goldhill.png'
      };
  end
  for ii = 1:length(listOfImages)
    iname = listOfImages{ii};
    I = imread(iname);
    [route,name,ext,ver] = fileparts(iname);
    if ~exist([route '/noisy'],'file')
      mkdir([route '/noisy']);
    end
    for si = 1:length(listOfSigmas)
      S = listOfSigmas(si);
      In = AddNoise(I,'Gaussian',S);
      iname = sprintf('%s/noisy/%s_g%02d%s',route,name,S,ext);
      imwrite(In,iname);
    end
  end
end
