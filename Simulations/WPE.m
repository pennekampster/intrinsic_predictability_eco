function [ H ] = WPE( x, mwl, tau )
%WPE: Garland weighted permutation entropy 
%   Estimates the entropy of any real-valued time series

H=zeros(1,mwl-2);

xl=x;
    for m=1:mwl-1
        xl=[xl vertcat(x(tau*m+1:end),zeros(tau*m,1))]; %table of lags
    end
    xl(end-(mwl*tau):end,:)=[];% remove rows with zeros at the end of xl ('end' is the last row, 'end'-1 is the row before...)    

    %loop over permutation lengths from 3 to maximum word length
       for wl=3:mwl 
            xltemp=xl(:,1:wl); %just the columns of xl for this word length
            maxH=sum(log2(2:wl)); %maximum entropy
            [~,ind]=sort(xltemp,2);
            word=ind*(10.^(0:wl-1)'); %list of words in the time series
            vword=var(xltemp,[],2); %calculate the weight of each word as the variance of the word elements
            u=unique(word);
            for k=1:length(u) %for each unique word
                inc=(word==u(k)); %index the matching words in the time series word list
                vw(k)=sum(inc.*vword); %sum the variances of that word
            end
          
            W=vw/sum(vw); %the weighted probability of a permutation
            
            %if any of the variances are zero then this returns W=NaN when it should be 1
            if vw==0
                W=1;
            end
                                   
            Hobs=nansum(W.*log2(W)); %the weighted permutation entropy with nansum incase one of the W is zero
            
            H(wl-2)=-Hobs/maxH; %normalized WPE relative to the maximum entropy
            clear vw k inc
       end

end

