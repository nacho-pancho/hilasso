function stego_dir=embed(cover_dir, stego_dir,img_fmt, embedder, rate, M)
system(["mkdir -p " stego_dir]);
printf("Stego images will be written in %s\n",stego_dir);      
printf("Creating %d stego images with embedding rate %f\n",M,rate);
images=filenames(cover_dir,img_fmt);
if nargin < 5
  M = size(images,1);
end

for i=1:M
  img_file=images(i).name;
  [fdir,fname,fext]=splitFileName(img_file);
  stego_file=[stego_dir '/' fname '.' fext];
  fstat=stat(stego_file);
  if isempty(fstat) % do not compute twice
    cover_file = [cover_dir '/' fname '.' fext];
    X=imread(cover_file);
    if embedder == "random"
      S=randembed(X,rate);
    else if embedder == "ternary"
      S=ternaryEmbed(X,rate);
    else
      S=mhpdmEmbed(X,rate);
    end
    imwrite(stego_file,S);
    continue;
  end
end
    