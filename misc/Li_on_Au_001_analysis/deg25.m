clear all
clc
close all

%{
I will work on multiple gaussian fits here
https://en.wikipedia.org/wiki/Full_width_at_half_maximum
%}
fid=fopen('Li_on_Au_001_25Deg_nofric.txt','r');


%the first line has the recorded names of values
names_of_values=textscan(fid, '%s %s %s %s ', 1);
%the second line has the data for these values
first_data_value=textscan(fid, '%f %f %f %f ',1);
remainingdata=textscan(fid, '%f %f ');

fclose (fid);

%now unpack the big stuff
energy(1,1)=first_data_value{1,1}(:,:);
%weight=every_data_value{1,8}(:,:);
intensity(1,1)=first_data_value{1,1}(:,:);

cnt0=1;

for jj = 2:1:1000
  energy(jj,1)=remainingdata{1,1}(cnt0,1);
  intensity(jj,1)=remainingdata{1,2}(cnt0,1);
  cnt0=cnt0+1;
 
endfor

%Now lets look for local max values
cnt=0;
localmax =0;
localindex=0;
localwidth=0;
for kk = 6:length(intensity)-6
  
  #i establish anything below 0.1 as the noise
  if intensity(kk,1) >= 0.1
    if   intensity(kk,1) > intensity(kk-1,1) && intensity(kk,1)> intensity(kk+1,1) 
    
      localmax(cnt+1,1)=intensity(kk,1)   ;
      localindex(cnt+1,1)=energy(kk,1);
      %update counter
      cnt=cnt+1;
    endif
  endif
  
endfor

cnt2=1; localvalsdown=0; cnt3=0; localindexdown=0;

%Now let's find linear values
 cnt6=1;
for mm = 4:length(intensity)-4
     %The problem wit this fwhm is due to the fact that the smaller peak does not actually have a datapoint at it's FWHM
     %This works for the larger gaussian peak
      for pp =1#:length(localmax)
        if ((intensity(mm,1) <= (localmax(pp)/2)) &&  (intensity(mm+1,1) >= (localmax(pp)/2))) || ((intensity(mm-1,1) >= (localmax(pp)/2)) &&  (intensity(mm,1) <= (localmax(pp)/2)))
          %we are approximately at the fullwidth half max location of peaks (the actual fwhm is every two data values in here)
          FWHM(cnt6,1)=energy(mm,1);
          cnt6=cnt6+1;
        endif
      endfor
      
    #i establish anything below 0.1 as the noise
  if intensity(mm,1) >= 0.1
    if  (intensity(mm-1,1)  < intensity(mm,1) ) &&  (intensity(mm,1) < intensity(mm+1,1)) 
      localval(cnt2,1) = intensity(mm,1);
      localindexup(cnt2,1) = energy(mm,1);
      %update counter
      cnt2=cnt2+1;

    elseif (intensity(mm,1) < intensity(mm-1,1)) && (intensity(mm,1)> intensity(mm+1,1)) 

     localvalsdown(cnt3+1 ,1 )=intensity(mm,1);
     localindexdown(cnt3+1,1)=energy(mm,1);
     cnt3=cnt3+1;
    
    endif
  endif
  
endfor
figure(1)
subplot(2,1,1)
hold on
plot(localindexup, localval , '.k')
plot(localindexdown, localvalsdown, '+r')
plot(localindex, localmax, 'm*')
xlabel('Energy eV')
ylabel('intensity')
subplot(2,1,2)
plot(energy,intensity, 'b.')
xlim([0.86 0.95])
ylabel('intensity')
xlabel('Energy eV')

%to get an approximate of the full width
cnt4=1;
for nn = 1:length(localvalsdown)-1
  %attempt to get the first peaks right most data point width
   if localvalsdown(nn,1) < localvalsdown(nn+1,1)
     %then we should be at the next set of values going update
     energysave(cnt4,1)=localindexdown(nn,1);
     cnt4=cnt4+1;
   endif
   
endfor

cnt5=1;
for oo = 1:length(localval)-1
  %attempt to get the second peak's leftmost datapoint
   if localval(oo,1) > localval(oo+1,1)
     %then we should be at the next set of values going update
     energysave1(cnt5,1)=localindexup(oo,1);
     cnt5=cnt5+1;
   endif
   
endfor

peak1width= 100*(FWHM(2,1)-FWHM(1,1));
peak2width= 100*(localindexdown(end,1)-energysave);

%recall FWHM = 2.35*sigma
    %https://en.wikipedia.org/wiki/Gaussian_function
sigma1=peak1width/2.35;
sigma2=peak2width/2.35;
amp1=localmax(1,1);
amp2=localmax(2,1);
x1=100*localindex(1,1)
x2=100*localindex(2,1);

fprintf('sigma 1 , sigma2, amp1, amp2, x1, x2 \n');
fprintf('%f %f %f %f %f %f \n', sigma1, sigma2, amp1, amp2, x1, x2);


