data_path = '';

traindata = cell(1,50);
trainX = cell(1,10);
forg_tr = cell(1,25);
gen_tr = forg_tr;

for i=0:0
        for j=0:0
            fname = sprintf('00%d%d',i,j);
            fname1 = [fname];
            dir1=dir(fname1);
            for k=3:52
                   filename = dir1(k).name;
                   filename1 = sprintf('00%d%d/%s',i,j,filename);
                   [x y z az in pps] = FPG_Signature_Read(filename1,0,0);
                   
                   si = size(x); %size of each vector
                   tstamp = zeros(si); %timestamp
                   for p=1:si
                        tstamp(p) = p/pps;
                   end
                   traindata{i*10+j+1,k-2} = [x y z az in tstamp];
                   if(k<28)
                       forg_tr{i*10+j+1,k-2}=[x y z az in tstamp];
                   else
                       gen_tr{i*10+j+1,k-27}=[x y z az in tstamp];
                   end
            end
        end
end


%disp(gen_tr{1,1});

% max no. of points needed for each signature after interpolation
max_gen_data_points = 0; 
% min no of points in genuine signatures of the given user.
min_gen_data_points = 99999;

for i=1:25
    [a b] = size(gen_tr{1,i});
    
    if(a > max_gen_data_points) 
        max_gen_data_points = a;
    end
    if(a < min_gen_data_points)
        min_gen_data_points = a;
    end
    
    % matrix storing values for a ith genuine signature
    gen_mat = gen_tr{1,i};  
    
    % calculating delta
    for j=2:a
       gen_mat(j-1,1:b-1) = gen_mat(j,1:b-1) - gen_mat(j-1,1:b-1);
       gen_mat(j-1,b) = gen_mat(j-1,b) - gen_mat(1,b);
    end
    % last row is redundant so removed
    gen_mat(a,:) = [];
    
    % updation of global matrix
    gen_tr{1,i} = gen_mat;
end

%disp(max_gen_data_points);
%disp(min_gen_data_points);

for i=1:25
    [a b] = size(forg_tr{1,i});
    
    % matrix storing values for ith forged signature
    forg_mat = forg_tr{1,i};
    
    % calculating delta for the forged signature
    for j=2:a
       forg_mat(j-1,1:b-1) = forg_mat(j,1:b-1) - forg_mat(j-1,1:b-1);
       forg_mat(j-1,b) = forg_mat(j-1,b) - forg_mat(1,b);
    end
    % last row is redundant so removed
    forg_mat(a,:) = [];
    
    % updating the global matrix for ith forged signature of user
    forg_tr{1,i} = forg_mat;
end

% For Normalization
% calculating mean and standard deviation for genuine signatures
for i=1:25
  gen_mat = gen_tr{1,i};
  a =  size(gen_mat);
  
  % mean_val contains mean for all features using all samples
  mean_val = zeros(5);
  for j=1:5
      mean_val(j) = mean(gen_mat(:,j));
  end
  
  % std_val contains std deviation for all features using all samples
  std_val = zeros(5);
  for j=1:5
      std_val(j) = std(gen_mat(:,j));
  end
  
  % undergoing normalization of ith genuine signature
  for j=1:5%columns(features)
     for k=1:a%rows
         gen_mat(k,j) = (gen_mat(k,j) - mean_val(j))/std_val(j);
     end
  end
  
  % updating the global matrix with normalized values
  gen_tr{1,i} = gen_mat;
end

% For Normalization
% calculating mean and standard deviation for forged signatures
for i=1:25
  forg_mat = forg_tr{1,i};
  a =  size(forg_mat);
  
  % mean_val contains mean for all features using all samples
  mean_val = zeros(5);
  for j=1:5
      mean_val(j) = mean(forg_mat(:,j));
  end
  
  % std_val contains std deviation for all features using all samples
  std_val = zeros(5);
  for j=1:5
      std_val(j) = std(forg_mat(:,j));
  end
  
  % undergoing normalization of ith forged signature
  for j=1:5%columns(features)
     for k=1:a%rows
         forg_mat(k,j) = (forg_mat(k,j) - mean_val(j))/std_val(j);
     end
  end
  
  % updating the global matrix with normalized values
  forg_tr{1,i} = forg_mat;
end

fprintf('Sizes of genuine signatures before interpolation');
for i=1:25
    disp(size(gen_tr{1,i}));
end


% For Interpolation
% interpolation of genuine signatures
for i=1:25
    gen_mat = gen_tr{1,i};
    
    %cur_size -> no of sample points presently in ith sign & k -> features
    [cur_size k] = size(gen_mat);
    
    t_stamp = gen_mat(:,k);
    
    % n -> extra points points to be added
    n = max_gen_data_points - cur_size;
    
    interval = floor(cur_size/(n+1));
    
    % new_tstamp --> column vector of n elements.
    extra_tstamp = zeros(n,1);
    
    for j=1:n
        extra_tstamp(j,1) = (t_stamp(j*interval,1) + t_stamp(j*interval + 1 ,1))/2;
    end
    
    % yi = interp1(x,Y,xi) Y is the function to be interpolated.
    % for each feature , need to get n new sample points.
    
    new_points = zeros(n,k);
    
    % for each feature j we are adding new_points.
    for j=1:k-1
        new_points(:,j) = interp1(t_stamp,gen_mat(:,j),extra_tstamp);
    end
    
    % appending extra_tstamp vector for sorting to be done later.
    new_points(:,k) = extra_tstamp;
    
    % concatenating of the new sample points to the original set of points
    new_gen_mat = [gen_mat ; new_points];
    
    % sorting the sample points based on timestamp values (k th col)
    new_sorted_gen_mat = sortrows(new_gen_mat,k);
    
    % updating the global array with values after interpolation
    gen_tr{1,i} = new_sorted_gen_mat;
end

fprintf('Sizes of genuine signatures after interpolation');
for i=1:25
    disp(size(gen_tr{1,i}));
end

% For Interpolation
% interpolation of forged signatures
for i=1:25
    forged_mat = forg_tr{1,i};
    
    %cur_size -> no of sample points presently in ith sign & k -> features
    [cur_size k] = size(forged_mat);
    
    t_stamp = forged_mat(:,k);
    
    % n -> extra points points to be added
    n = max_gen_data_points - cur_size;
    
    % if n < 0 remove unwanted rows from forged_mat
    if n < 0
       for j = cur_size : -1 : max_gen_data_points + 1
           %disp(n); disp(forged_mat(cur_size,:));
           forged_mat(j,:) = [];
       end   
    elseif n ~= 0
        interval = floor(cur_size/(n+1));
    
        % new_tstamp --> column vector of n elements.
        extra_tstamp = zeros(n,1);
    
        for j=1:n
            extra_tstamp(j,1) = (t_stamp(j*interval,1) + t_stamp(j*interval + 1 ,1))/2;
        end
    
        % yi = interp1(x,Y,xi) Y is the function to be interpolated.
        % for each feature , need to get n new sample points.
    
        new_points = zeros(n,k);
    
        % for each feature j we are adding new_points.
        for j=1:k-1
            new_points(:,j) = interp1(t_stamp,forged_mat(:,j),extra_tstamp);
        end
    
        % appending extra_tstamp vector for sorting to be done later.
        new_points(:,k) = extra_tstamp;
    
        % concatenation of the new sample points to the original set of points
        new_forged_mat = [forged_mat ; new_points];
    
        % sorting the sample points based on timestamp values (k th col)
        new_sorted_forged_mat = sortrows(new_forged_mat,k);
        
        forged_mat = new_sorted_forged_mat;
    end
        % updating the global array with values after interpolation
        forg_tr{1,i} = forged_mat;
end



fprintf('Sizes of forged signatures after interpolation');
for i=1:25
    disp(size(forg_tr{1,i}));
end

for i=1:5
    trainX{1,i} = gen_tr{1,i};
end

for i=1:5
    trainX{1,5+i} = forg_tr{1,i};
end

trainY = [1 1 1 1 1 -1 -1 -1 -1 -1];

%calculation of threshold(e)

g_max = 0;
max = 0;
for i=1:5
    for j=i+1:5
        min = realmax;
        max = 0;
        X = trainX{1,i};
        Y = trainX{1,j};
        [si_y k] = size(Y);
        [si_x k] = size(X);
        for k=1:si_x
            for l=1:si_y
                m = X(k,:);
                n = Y(l,:);
                d = dist(m,n);
                if(d<min)
                    min = d;
                end
            end
            if(min > max)
                max = min;
            end
        end
    end
    if(g_max < max)
        g_max = max;
    end
end

thres = g_max;

%disp(thres);

%the kernel

%md1 = fitcsvm(trainX,trainY,'KernelFunction','edrkernel','Standardize',true);