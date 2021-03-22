function images = load_images(path,wildcard)
% images = load_images(path,wildcard) returns the a structure with the 
% image data of all images in the defined path that match the wildcard

img_files = dir([path wildcard]);

for i = 1:length(img_files)
    images{i} = imread([path img_files(i).name]);
end