function [g, gd, regmax, err] = gettf1(u, y, nn, tt, flag)

%[g, gd, regmax] = gettf1(u,y,nn,tt,flag)
%
% 'gettf1' obtains a transfer function
% from a time sequence and regressor
% g --> Estimated continuous transfer function Output/Input (y/u)
% gd--> Estimated discrete transfer function Output/Input (y/u)
% regmax --> Best regressor when used with optimization (i.e. flag=1)
% u --> Input sequence
% y --> Output sequence
% nn --> Regressor [na nb nk]
% tt --> Time sequence
% flag --> (0 to just use the given regressor, and 1 to search for
%           the optimum regressor using the given one as the upper bound)
% Author: Y. Zeyada 16 Aug 2000

% The next 7 lines were added by Ron Hess 11/9/01
ts=tt(6)-tt(5);
[bin,ain]=butter(4,6*ts);
[bout,aout]=butter(4,6*ts);
in=filter(bin,ain,u);
out=filter(bout,aout,y);
u=in;
y=out;

u=u';
y=y';


NN=nn;

if flag==1,
    val=[];
    regr=[];
    for NA=1:nn(1),
        for NB=1:nn(2),
            for NK=-3:nn(3),
                na=NA;
                nb=NB;
                nk=NK;
                j=0;
                a=[];
                b=[];
                for i=10+(nk-1):length(u)-(na+(nk-1))-10,
                    j=j+1;
                    a(j,1:na)=y(i:i+na-1);
                    a(j,na+1:na+nb)=u(i+na-nk-nb+1:i+na-nk);
                    b(j,1)=y(i+na);
                end
                x=a\b;
                xy=x(na:-1:1);
                xu=x(na+nb:-1:na+1);
                xy1=zeros(1,nk-1);

                gdd=tf([xu'],[1 -xy' xy1 ],ts);
                gg=d2c(gdd,'tustin');
                gg=minreal(gg);

                [yc,t,x]=lsim(gg,u,tt);

                err=yc'-y;
                ems=sqrt(err*err');
                val=[val;1/(ems+eps)];
                regr=[regr;NA NB NK];
            end
        end
    end
    maxval=max(val);
    imaxval=find(val==maxval);
    NN=regr(imaxval,:);
end

a=[];
b=[];

if isempty(NN)==0,
    na=NN(1);
    nb=NN(2);
    nk=NN(3);
    j=0;
    for i=10+(nk-1):length(u)-(na+(nk-1))-10,
        j=j+1;
        a(j,1:na)=y(i:i+na-1);
        a(j,na+1:na+nb)=u(i+na-nk-nb+1:i+na-nk);
        b(j,1)=y(i+na);
    end
    x=a\b;
    xy=x(na:-1:1);
    xu=x(na+nb:-1:na+1);
    xy1=zeros(1,nk-1);

    gdd=tf([xu'],[1 -xy' xy1 ],ts);
    gg=d2c(gdd,'tustin');
    gg=minreal(gg);

    [yc,t,x]=lsim(gg,u,tt);
    plot(tt,y,tt,yc','r')
    pause;
    err=yc'-y;
    plot(tt,err)
    rmsy=sqrt(y*y'*ts/max(tt));
    rmserr=sqrt(err*err'*ts/max(tt))/rmsy
    ems=sqrt(err*err');
    val=1/(ems+eps)
    gd=gdd;
    g=gg;
    regmax=NN;
else
    g=[];
    gd=[];
    regmax=[];
end
