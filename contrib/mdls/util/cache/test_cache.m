more off

fdir='/home/nacho/images';
fext='png';
fcache='numpix';

images=filenames(fdir,fext)
numpixels=zeros(length(images),1);
N=length(images);

cache_data=cache_load(fdir,fext,fcache);

for i=1:1:N% exclude . and ..
  img_name = images(i).name;
  entry=cache_read(fdir,fext,cache_data,img_name);
  printf("%s",img_name);
  if (isempty(entry))
    X=imread([fdir '/' img_name]);
    numpixels(i)=size(X,1)*size(X,2);
    cache_data=cache_write(fdir, fext, cache_data, img_name, numpixels(i));
    printf(":%d [computed]\n",numpixels(i));
  else
    numpixels(i)=entry;
    printf(":%d [cached]\n",numpixels(i));
  end
end

for i=3:length(images)
  img_name = images(i).name;
  entry=cache_read(fdir,fext,cache_data,img_name);
  printf("%s",img_name);
  if (isempty(entry))
    X=imread([fdir '/' img_name]);
    numpixels(i)=size(X,1)*size(X,2);
    cache_data=cache_write(fdir,fext,cache_data,img_name,numpixels(i));
    printf(":%d [computed]\n",numpixels(i));
  else
    numpixels(i)=entry;
    printf(":%d [cached]\n",numpixels(i));
  end
end
cache_save(fdir,fext,cache_data,fcache);
% status
cache_status(fdir,fext,cache_data)
