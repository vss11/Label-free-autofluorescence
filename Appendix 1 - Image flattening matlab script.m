clc
clear
close all
tic
currentFolder=pwd;

[FileName,PathName] = uigetfile('*.hdf5;*.h5','Select the file');
cd(currentFolder)
cd(PathName)
%
trial=0;
optimize=0;


info=h5info(strcat(PathName,FileName));
dataset=info.Datasets.Name;
rawdata=h5read(strcat(PathName,FileName),['/' dataset]);            %%% original
%rawdata(:,:,:,1)=h5read(strcat(PathName,FileName),['/' dataset{1}]);
%rawdata(:,:,:,2)=h5read(strcat(PathName,FileName),['/' dataset{2}]);
%rawdata(:,:,:,3)=h5read(strcat(PathName,FileName),['/' dataset{3}]);


%% select shg
timepoint=size(rawdata,5);
allheight=cell(timepoint,1);
rawdata=rawdata(1,:,:,:,:);                     %%% SHG is channel 1
rawdata=permute(squeeze(rawdata),[2 1 3 4]);    %%% ? why [2,1,3,4]

%%
close all
trial=0;
for t=1:timepoint
    while(1)
% gather initial parameters        
        prompt = {'Enter XY resolution (um/pixel):','Enter Z resolution (um/pixel):',...
            'Enter approx. frame # where BM starts','XY blurring (in um, typically 10-25)',...
            'Z blurring (in um, typically 1-4)','Downsample resolution (um/pixel)',...
            '# of frames in each zstack','# of timepoints per z stack'};
        dlg_title = 'Input';
        num_lines = 1;
        
        
        
        if trial ==0 && t==1
    
                defaultans = {'0.3','0.4',...
                    '60','20',...
                    '1.5','1',...
                    '74',num2str(timepoint)};                  

            
        else
            defaultans = {num2str(res_xy), num2str(res_z),...
                num2str(BM_start), num2str(blur),...
                num2str(blurz), num2str(dsrate),...
                num2str(stacksize),num2str(timepoint)};
        end
        
        
        if optimize == 0
            answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
            res_xy=str2num(answer{1}); %XY RESOLUTION
            res_z=str2num(answer{2}); %Z RESOLUTION
            BM_start=str2num(answer{3}); %ROUGHLY WHERE THE BASEMENT MEMBRANE STARTS
            blur=str2num(answer{4}); %AMOUNT TO BLUR
            blurz=str2num(answer{5}); %AMOUNT TO BLUR
            dsrate=str2num(answer{6}); %Down sampled resolution
            stacksize=str2num(answer{7}); %number of slices in the zstack
            timepoint=str2num(answer{8}); %SHG channel number
            
        else if optimize == 1
            end
        end
          
        
        
% extract out SHG 
        blue=squeeze(rawdata(:,:,:,t));              
        stacksize=size(blue,3);
        movie=zeros(size(blue,1),size(blue,2),3,stacksize);    %%%% why need a 4D matrix?
        
        movie(:,:,1,:)=blue;
        movie(:,:,2,:)=blue;
        movie(:,:,3,:)=blue;

        
        originalsize_x=size(movie,1);
        originalsize_y=size(movie,2);
        dsrate_adjust=(res_xy/dsrate);
        if originalsize_x*dsrate_adjust > floor(originalsize_x*dsrate_adjust)
        dssize_x=floor(originalsize_x*dsrate_adjust)+1;
        else
         dssize_x=floor(originalsize_x*dsrate_adjust);
        end
        if originalsize_y*dsrate_adjust > floor(originalsize_y*dsrate_adjust)
        dssize_y=floor(originalsize_y*dsrate_adjust)+1;
        else
         dssize_y=floor(originalsize_y*dsrate_adjust);  
        end
        
        zslices=size(movie,4);
        filtSHG=zeros(dssize_x,dssize_y,zslices);
        
        
% call a function to first get rid of the cell signal from SHG
        
         movie(:,:,3,:) = SHGclean(squeeze(movie(:,:,3,:)),300);
        
        
%% blur the SHG

        movie(:,:,3,1:BM_start) = 0;
        
        for z=1:zslices %%% can consider do z= BMStart:zslice
            filtSHG(:,:,z)=imgaussfilt(imresize(movie(:,:,3,z),dsrate_adjust),blur);        %%% chose to blur the 3rd channel
            
            j=figure(1); clf;
            set(gcf,'Position',[225 420 1015 415]);        
            subplot(1,2,1);
            titlelabel=sprintf('Timepoint %d / %d; plane %d / %d',t,timepoint,z,zslices);
            imshow(filtSHG(:,:,z),[],...
                'InitialMagnification','fit');
            title(titlelabel);
            subplot(1,2,2);
            imshow(blue(:,:,z),[],'InitialMagnification','fit');
            title(titlelabel);
            
            pause(0.01);
        end
        
    
        threedblur=blurz/res_z;
        filtSHG=imgaussfilt3(filtSHG,threedblur);
        height=zeros(dssize_x,dssize_y);
        
% find SHG peaks in the image to generate height map        
        parfor x=1:dssize_x
            for y=1:dssize_y
                locs=[];
                test=squeeze(filtSHG(x,y,BM_start:end));
                % [peak,locs]=findpeaks(test,'MinPeakProminence',3,'SortStr','descend');
                [peak,locs]=findpeaks(test);
                [mtest,BM_end]=max(test);
                %  [~,locs]=findpeaks(test);
                if isempty(locs) == 1
                    if (BM_end + BM_start) > size(filtSHG,3)-10
                        height(x,y)=(BM_end + BM_start);
                    else
                        height(x,y)=BM_start;
                        continue
                    end
                else
                    if mtest > 1.5*peak
                        height(x,y)=(BM_end + BM_start);
                    else
                        [~,mp]=max(peak);
                        height(x,y)=locs(mp)+BM_start;
                    end
                end
            end
        end
        
% remove holes, smooth the height map       
        
        boxsize=round(dssize_x/5);
        height2=height;
        %remove outliers
        for x=1:dssize_x
            for y=1:dssize_y
                if height(x,y)==BM_start
                    %check limits
                    if x<(boxsize+1)
                        minx=1;
                    else
                        minx=x-boxsize;
                    end
                    
                    if x>(dssize_x-boxsize)
                        maxx=dssize_x;
                    else
                        maxx=x+boxsize;
                    end
                    
                    if y<(boxsize+1)
                        miny=1;
                    else
                        miny=y-boxsize;
                    end
                    
                    if y>(dssize_y-boxsize)
                        maxy=dssize_y;
                    else
                        maxy=y+boxsize;
                    end
                    
                    height2(x,y)=median(median(height(minx:maxx,miny:maxy)));
                    
                end
            end
        end
        
% filter and display height map     
        h=figure(2);
        height3=medfilt2(height2,[blur blur],'symmetric');
        height4=imresize(height3,1/dsrate_adjust);
        height5=imgaussfilt(height4,blur/dsrate_adjust);
        height6=round(height5);
        imagesc(height6);
        c=colorbar;
        ylabel(c,'z-slice');
        axis off
        titlelabel2=sprintf('Timepoint %d / %d',t,timepoint);
        title(titlelabel2);
        
        allheight{t}=height6;
        
        if optimize ==0
            promptMessage = sprintf('Do you want to continue with this height map ,\n or try with different parameters?');
            button = questdlg(promptMessage, 'Continue?', 'Continue', 'Restart', 'Continue');
        end
        
        if strcmpi(button, 'Continue')
            optimize=1;
            break
        end
        trial=trial+1;
    end

end



%% plot curvature of z stack

%figure(3); clf;

for i=1:timepoint
    titlelabel=sprintf('Timepoint %d / %d',i,timepoint);
    figure(1); clf;
    min_height=min(min(allheight{i}));
    max_height=max(max(allheight{i}));
    h=surf(allheight{i},'FaceLighting','gouraud','FaceColor','interp');
    set(h,'LineStyle','none'); 
    view([-19,15]);
    zlim([30 55]);
    title(titlelabel);
    pause(0.01);
% figure(2);
% plot(i,max_height-min_height,'o','MarkerSize',5); hold on
% xlabel('Time');
% ylabel('Max height - min height');
end


%% IF HEIGHTMAP IS GOOD, SAVE allheight as allheight.mat %%%%%%%%%%%%
save('allheight.mat');
%%
%NORMALIZE AND SAVE DATA 

close all

load allheight.mat
data=[];
data=h5read(FileName,'/data');              %% original data=h5read('FinalAnalysis.hdf5','/data');
%data=h5read(FileName,'/data');
%data(:,:,:,1)=h5read(strcat(PathName,FileName),['/' dataset{1}]);
%data(:,:,:,2)=h5read(strcat(PathName,FileName),['/' dataset{2}]);
%data(:,:,:,3)=h5read(strcat(PathName,FileName),['/' dataset{3}]);

data=permute(data,[1 3 2 4 5]);  %%% original
%data=permute(data,[3 1 2 4 5]);     %%% updated on 06/08//2021
  
%%
heightstack=zeros(size(data,2),size(data,3),size(data,5));
for i=1:length(allheight)
    heightstack(:,:,i)=imresize(allheight{i},[size(data,2) size(data,3)]);
end

%% normalize
xsize=size(data,2);
ysize=size(data,3);
Num_channel = size(data,1);

for t=1:timepoint
    movie=squeeze(data(:,:,:,:,t));
    movie=permute(movie,[2 3 1 4]);    %%% original
    %movie=permute(movie,[2 3 4 1]);     %%% update 20210608
    height6=round(heightstack(:,:,t));
    %%% limit max height to not exceed stack
    height6(height6>stacksize) = stacksize;
    
    
    maxheight=max(max(height6));
    minheight=min(min(height6));
    heightdiff=maxheight-minheight;
    %%newzslices=70+heightdiff;                   %%%%% original
    newzslices=stacksize+heightdiff;
    newstack=zeros(xsize,ysize,Num_channel,newzslices);   %%%% the third dimension should be the number of channesl
    
    %pad zeros
    movie2=cat(4,zeros(xsize,ysize,Num_channel,heightdiff),movie);    %% the third numebr should be the number of channals
    movie3=cat(4,movie2,zeros(xsize,ysize,Num_channel,heightdiff));
    
    f=waitbar(1,'Normalizing z-stack...');
    
    d = size(movie3,4);
    for x=1:xsize
        waitbar(x/xsize);
        for y=1:ysize
            newstack(x,y,:,1:(stacksize+maxheight-height6(x,y)))=movie3(x,y,:,(height6(x,y)-minheight+1):(heightdiff+stacksize));
            %%%newstack(x,y,:,1:(d-(height6(x,y)-minheight)))=movie3(x,y,:,heightdiff:(heightdiff+stacksize));                    
            %%%newstack(x,y,:,:)=movie3(x,y,:,height6(x,y)-minheight+[1:newzslices]);
        end
    end
    
    close(f);
    
    
    for z=1:newzslices
        zprofile(z)=sum(sum(sum(newstack(:,:,2,z))));           %%%% shouldn't this be only on SHG channal?
    end
    [locmax,loc]=max(zprofile);
    
    locmin=5; %how many um below the zero plane to start (i.e. basement membrane)
    locmax=loc+find(zprofile(loc:end)<0.25*locmax,1);
    if isempty(locmax)
        locmax = newzslices;
    end
    allmovies{t}=newstack(:,:,:,loc-round(locmin/res_z):locmax-1);
    
    %%%% don't know why so skip this
    
    allmovies{t}=permute(newstack,[3 1 2 4]);
    
    [a,b,c,d]=size(allmovies{t});
    
    if t >1
        while size(allmovies{t},4) > size(allmovies{1},4);
            allmovies{t}(:,:,:,end)=[];
        end
        while size(allmovies{t},4) < size(allmovies{1},4);
            allmovies{t}(:,:,:,end+1)=zeros(a,b,c);
        end
    end
    
    if t==1
        allmoviestiff=zeros([size(allmovies{1}) timepoint],'uint16');
        allmoviestiff(:,:,:,:,t)=uint8(allmovies{t});
        %allmoviestiff(:,:,:,:,t)=allmovies{t};
    end
    
    allmoviestiff(:,:,:,:,t)=uint16(allmovies{t});
    %allmoviestiff(:,:,:,:,t)=allmovies{t};
    t
end

%%


savename='normStain.h5';
h5create(savename,'/data',[size(allmoviestiff)],'Datatype','uint16');
h5write(savename,'/data',uint16(allmoviestiff));
%h5create(savename,'/data',[size(allmoviestiff)]);
%h5write(savename,'/data',allmoviestiff);
