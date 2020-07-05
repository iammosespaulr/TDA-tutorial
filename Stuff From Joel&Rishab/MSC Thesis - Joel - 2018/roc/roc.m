function ROCout=roc(x,varargin)
% ROC - Receiver Operating Characteristics.
% The ROC graphs are a useful tecnique for organizing classifiers and
% visualizing their performance. ROC graphs are commonly used in medical
% decision making.
%
% Syntax: ROCout=roc(x,thresholds,alpha,verbose)
%
% Input: x - This is a Nx2 data matrix. The first column is the column of the data value;
%            The second column is the column of the tag: unhealthy (1) and
%            healthy (0).
%        Thresholds - If you want to use all unique values in x(:,1)
%            then set this variable to 0 or leave it empty;
%            else set how many unique values you want to use (min=3);
%        alpha - significance level (default 0.05)
%        verbose - if you want to see all reports and plots (0-no; 1-yes by
%        default);
%
% Output: if verbose = 1
%         the ROCplots, the sensitivity and specificity at thresholds; the Area
%         under the curve with Standard error and Confidence interval and
%         comment.
%         if ROCout is declared, you will have a struct:
%         ROCout.AUC=Area under the curve (AUC);
%         ROCout.SE=Standard error of the area;
%         ROCout.ci=Confidence interval of the AUC
%         ROCout.co=Cut off points
%         ROCdata.xr and ROCdata.yr points for ROC plot
%
%           Created by Giuseppe Cardillo
%           giuseppe.cardillo-edta@poste.it
%
% To cite this file, this would be an appropriate format:
% Cardillo G. (2008) ROC curve: compute a Receiver Operating Characteristics curve.
% http://www.mathworks.com/matlabcentral/fileexchange/19950

%Input Error handling
p=inputParser;
addRequired(p,'x',@(x) validateattributes(x,{'numeric'},{'2d','real','finite','nonnan','nonempty','ncols',2}));
addOptional(p,'threshold',0, @(x) isnumeric(x) && isreal(x) && isfinite(x) && isscalar(x) && (x==0 || x>2));
addOptional(p,'alpha',0.05, @(x) validateattributes(x,{'numeric'},{'scalar','real','finite','nonnan','>',0,'<',1}));
addOptional(p,'verbose',1, @(x) isnumeric(x) && isreal(x) && isfinite(x) && isscalar(x) && (x==0 || x==1));
parse(p,x,varargin{:});
threshold=p.Results.threshold; alpha=p.Results.alpha; verbose=p.Results.verbose;
clear p
assert(all(x(:,2)==0 | x(:,2)==1),'Warning: all x(:,2) values must be 0 or 1')
if all(x(:,2)==0)
    error('Warning: there are only healthy subjects!')
end
if all(x(:,2)==1)
    error('Warning: there are only unhealthy subjects!')
end

tr=repmat('-',1,100);
z=sortrows(x,1);
if threshold==0
    labels=unique(z(:,1));%find unique values in z
else
    K=linspace(0,1,threshold+1); K(1)=[];
    labels=quantile(unique(z(:,1)),K)';
end
clear z

labels(end+1)=labels(end)+1;
ll=length(labels); %count unique value
a=zeros(ll,2); c=zeros(ll,1);%array preallocation
ubar=mean(x(x(:,2)==1),1); %unhealthy mean value
hbar=mean(x(x(:,2)==0),1); %healthy mean value
N=length(x);
for K=1:ll
    if hbar<ubar
        TP=length(x(x(:,2)==1 & x(:,1)>labels(K)));
        FP=length(x(x(:,2)==0 & x(:,1)>labels(K)));
        FN=length(x(x(:,2)==1 & x(:,1)<=labels(K)));
        TN=length(x(x(:,2)==0 & x(:,1)<=labels(K)));
    else
        TP=length(x(x(:,2)==1 & x(:,1)<labels(K)));
        FP=length(x(x(:,2)==0 & x(:,1)<labels(K)));
        FN=length(x(x(:,2)==1 & x(:,1)>=labels(K)));
        TN=length(x(x(:,2)==0 & x(:,1)>=labels(K)));
    end
    M=[TP FP;FN TN];
    a(K,:)=diag(M)'./sum(M); %Sensitivity and Specificity
    c(K)=trace(M)/N; %Efficiency
end
clear ll K TP FP FN TN M N
b=[a(:,1)./(1-a(:,2)) (1-a(:,1))./a(:,2)];

xroc=1-a(:,2); yroc=a(:,1); %ROC points
if hbar>ubar
    xroc=flipud(xroc); yroc=flipud(yroc); %ROC points
end
clear ubar hbar

st=[1 mean(xroc) 1]; L=[0 0 0]; U=[Inf 1 Inf];
fo_ = fitoptions('method','NonlinearLeastSquares','Lower',L,'Upper',U,'Startpoint',st);
ft_ = fittype('1-1/((1+(x/C)^B)^E)','dependent',{'y'},'independent',{'x'},'coefficients',{'B', 'C', 'E'});
rocfit = fit(xroc,yroc,ft_,fo_);


xfit=linspace(0,1,500);
yfit=rocfit(xfit);
clear st L U fo_ ft_ 
Area=trapz(xfit,yfit); %estimate the area under the curve
%standard error of area
lu=length(x(x(:,2)==1)); %number of unhealthy subjects
lh=length(x(x(:,2)==0)); %number of healthy subjects
Area2=Area^2; Q1=Area/(2-Area); Q2=2*Area2/(1+Area);
V=(Area*(1-Area)+(lu-1)*(Q1-Area2)+(lh-1)*(Q2-Area2))/(lu*lh);
Serror=realsqrt(V);
%confidence interval
ci=Area+[-1 1].*(realsqrt(2)*erfcinv(alpha)*Serror);
if ci(1)<0; ci(1)=0; end
if ci(2)>1; ci(2)=1; end
%z-test
SAUC=(Area-0.5)/Serror; %standardized area
p=1-0.5*erfc(-SAUC/realsqrt(2)); %p-value

if nargout
    ROCout.AUC=Area; %Area under the curve
    ROCout.SE=Serror; %standard error of the area
    ROCout.ci=ci; % 95% Confidence interval
    ROCout.xr=xroc; %graphic x points
    ROCout.yr=yroc; %graphic y points
    ROCout.rocfit=rocfit;
end

clear lu lh Area2 Q1 Q2 V rocfit

if verbose==1
    %Performance of the classifier
    if Area==1
        str='Perfect test';
    elseif Area>=0.90 && Area<1
        str='Excellent test';
    elseif Area>=0.80 && Area<0.90
        str='Good test';
    elseif Area>=0.70 && Area<0.80
        str='Fair test';
    elseif Area>=0.60 && Area<0.70
        str='Poor test';
    elseif Area>=0.50 && Area<0.60
        str='Fail test';
    else
        str='Failed test - less than chance';
    end
    
    %display results
    disp('ROC CURVE ANALYSIS')
    disp(tr)
    disp(cell2table({Area, Serror, ci, str, SAUC, p},'VariableNames',{'AUC','Standard_error','Confidence_interval','Comment','Standard_AUC','p_value'}))
    if p<=alpha
        disp('The area is statistically greater than 0.5')
    else
        disp('The area is not statistically greater than 0.5')
    end
    disp(' ')
    %display graph
    H=figure;
    set(H,'Position',[4 402 560 420])
    axis square; hold on
    plot(xfit,yfit,'marker','none','linestyle','-','color','k','linewidth',2);
    patch(xfit,yfit,'g','FaceAlpha',0.5)
    patch([0 1 1],[0 0 1],'r','FaceAlpha',0.5)
    set(gca,'Xtick',0:0.1:1)
    grid on
    plot(xfit,yfit,'marker','none','linestyle','-','color','r','linewidth',2);
    H1=plot(xroc,yroc,'bo');
    set(H1,'markersize',6,'markeredgecolor','b','markerfacecolor','b')
    hold off
    xlabel('False positive rate (1-Specificity)')
    ylabel('True positive rate (Sensitivity)')
    title(sprintf('ROC curve (AUC=%0.4f)',Area))

    clear Area Serror ci str SAUC Serror xroc yroc H H1 xfit yfit
    if p<=alpha
        clear p
        ButtonName = questdlg('Do you want to input the true prevalence?', 'Prevalence Question', 'Yes', 'No', 'No');
        if strcmp(ButtonName,'Yes')  
            ButtonName = questdlg('Do you want to input the true prevalence as:', 'Prevalence Question', 'Ratio', 'Probability', 'Ratio');
            switch ButtonName
                case 'Ratio'
                    prompt={'Enter the Numerator or the prevalence ratio:','Enter the denominator or the prevalence ratio:'};
                    name='Input for Ratio prevalence';
                    Ratio=str2double(inputdlg(prompt,name));
                    POD=Ratio(1)/diff(Ratio); %prior odds
                case 'Probability'
                    prompt={'Enter the prevalence probability comprised between 0 and 1:'};
                    name='Input for prevalence';
                    pr=str2double(inputdlg(prompt,name));
                    POD=pr/(1-pr); %prior odds
            end
            clear ButtonName prompt name Ratio pr
            d=[1./(1+1./(b(:,1).*POD)) 1./(1+(b(:,2).*POD))];
            d((a(:,1)==0 & a(:,2)==1),1)=NaN;
            d((a(:,1)==1 & a(:,2)==0),2)=NaN;
            matrix=[labels'; a(:,1)'; a(:,2)';c';b(:,1)'; b(:,2)';d(:,1)'.*100; d(:,2)'.*100;]';
            clear POD a b c d labels
            if verbose
                disp('ROC CURVE DATA')
                disp(tr)
                disp(array2table(matrix,'VariableNames',{'Cut_off','Sensitivity','Specificity','Efficiency','PLR','NLR','Pos_pred','Neg_Pred'}))
            end
        else
            clear ButtonName
            matrix=[labels'; a(:,1)'; a(:,2)';c';b(:,1)'; b(:,2)']';
            clear a b c d labels
            if verbose
                disp('ROC CURVE DATA')
                disp(tr)
                disp(array2table(matrix,'VariableNames',{'Cut_off','Sensitivity','Specificity','Efficiency','PLR','NLR'}))
            end
        end
        clear tr
        
        CSe=matrix(matrix(:,1)==min(matrix(matrix(:,2)==max(matrix(:,2)))),1); %Max sensitivity cut-off
        CSp=matrix(matrix(:,1)==max(matrix(matrix(:,3)==max(matrix(:,3)))),1); %Max specificity cut-off
        CEff=matrix(matrix(:,1)==min(matrix(matrix(:,4)==max(matrix(:,4)))),1); %Max efficiency cut-off
        CPlr=matrix(matrix(:,5)==max(matrix(~isnan(matrix(:,5)) & ~isinf(matrix(:,5)),5)),1); %Max PLR cut-off
        CNlr=matrix(matrix(:,6)==max(matrix(~isnan(matrix(:,6)),6)),1); %Max NLR cut-off
        
        z=min(matrix(:,1));
        if z<0
            COEFF=abs(z)+1;
        else
            COEFF=0;
        end
        clear z
        
        M=mean(matrix(:,1))+COEFF; MM=max(matrix(:,1))+COEFF;
        ft_ = fittype('1-1/((1+(x/C)^B)^E)','dependent',{'y'},'independent',{'x'},'coefficients',{'B', 'C', 'E'});
        st=[1 M 1]; L=[0 0 0]; U=[Inf MM Inf];
        fo_ = fitoptions('method','NonlinearLeastSquares','Lower',L,'Upper',U,'Startpoint',st);
        fitSe = fit(matrix(:,1)+COEFF,matrix(:,2),ft_,fo_);
        
        st=[-1 M 1]; L=[-Inf 0 0]; U=[0 MM Inf];
        fo_ = fitoptions('method','NonlinearLeastSquares','Lower',L,'Upper',U,'Startpoint',st);
        fitSp=fit(matrix(:,1)+COEFF,matrix(:,3),ft_,fo_);
        
        ft_ = fittype('D+(A-D)/((1+(x/C)^B)^E)','dependent',{'y'},'independent',{'x'},'coefficients',{'A', 'B', 'C', 'D', 'E'});
        st=[min(matrix(:,4)) 1 M max(matrix(:,4)) 1]; L=[0 0 0 0 0]; U=[1 Inf MM 1 Inf];
        fo_ = fitoptions('method','NonlinearLeastSquares','Lower',L,'Upper',U,'Startpoint',st);
        fitEff=fit(matrix(:,1)+COEFF,matrix(:,4),ft_,fo_);

        myfun=@(x,se,sp) (-1./(1+(x./se(2)).^se(1)).^se(3))+(1./(1+(x./sp(2)).^sp(1)).^sp(3));
        SeSp=fzero(@(x) myfun(x,coeffvalues(fitSe),coeffvalues(fitSp)),M);
        
        clear M MM ft_ st fo_ myfun L U st

        xg=linspace(0,max(matrix(:,1))+COEFF,500);
        H2=figure;
        set(H2,'Position',[570 402 868 420])
        hold on
        H=ones(1,9);
        c=[0 0 1;1 0 0; 0 1 0; 0 0 0.1724; 1 0.1034 0.7241; 1 0.8276 0];
        H(1) = plot(xg,feval(fitSe,xg),'marker','none','linestyle','-','color',c(2,:),'linewidth',2);
        H(2)=plot([CSe CSe]+COEFF,[0 1],'marker','none','linestyle','--','color',c(2,:),'linewidth',2);
        H(3) = plot(xg,feval(fitSp,xg),'marker','none','linestyle','-','color',c(3,:),'linewidth',2);
        H(4)=plot([CSp CSp]+COEFF,[0 1],'marker','none','linestyle','--','color',c(3,:),'linewidth',2);
        H(5) = plot(xg,feval(fitEff,xg),'marker','none','linestyle','-','color',c(1,:),'linewidth',2);
        H(6)=plot([CEff CEff]+COEFF,[0 1],'marker','none','linestyle','--','color',c(1,:),'linewidth',2);
        H(7)=plot([SeSp SeSp],[0 1],'marker','none','linestyle','--','color',c(6,:),'linewidth',2);
        H(8)=plot([CPlr CPlr]+COEFF,[0 1],'marker','none','linestyle','--','color',c(5,:),'linewidth',2);
        H(9)=plot([CNlr CNlr]+COEFF,[0 1],'marker','none','linestyle','--','color',c(4,:),'linewidth',2);
        xlabel('Test cut-off')
        ylabel('Percent')
        hold off
        legend(H,...
            'Sensitivity',sprintf('Max Sensitivity cutoff: %0.4f',CSe),...
            'Specificity',sprintf('Max Specificity cutoff: %0.4f',CSp),...
            'Efficiency',sprintf('Max Efficiency cutoff: %0.4f',CEff),...
            sprintf('Cost Effective cutoff: %0.4f',SeSp-COEFF),...
            sprintf('Max PLR: %0.4f',CPlr),sprintf('Max NLR: %0.4f',CNlr),...
            'Location','BestOutside')
        axis([xg(1) xg(end) 0 1.1])
        if COEFF~=0
            xt=get(gca,'XTick'); Lxt=length(xt);
            xtl=cell(1,Lxt);
            for I=1:Lxt
                xtl{I}=sprintf('%0.2f',xt(I)-COEFF);
            end
            set(gca,'XTick',xt,'XTickLabel',xtl)
            clear xt xtl I Lxt
        end
        
        fprintf('1) Max Sensitivity Cut-off point= %0.4f\n',CSe)
        fprintf('2) Max Specificity Cut-off point= %0.4f\n',CSp)
        fprintf('3) Cost effective Cut-off point (Sensitivity=Specificity)= %0.4f\n',SeSp-COEFF),
        fprintf('4) Max Efficiency Cut-off point= %0.4f\n',CEff)
        fprintf('5) Max PLR Cut-off point= %0.4f\n',CPlr)
        fprintf('6) Max NLR Cut-off point= %0.4f\n',CNlr)
        m=[CSe CSp SeSp CEff CPlr CNlr];
        
    else
        m=NaN;
        matrix=NaN;
    end
    
    if nargout
        ROCout.co=m; % cut off points
        ROCout.table=matrix;
    end
end