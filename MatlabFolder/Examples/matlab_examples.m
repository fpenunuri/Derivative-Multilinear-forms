addpath('../01DualNumbersF;../02DMLF')

%example of section 5.1
m=20;
x=1:m;
y=sin(x);
z=cos(x);
w=sqrt(x);
q=log(x);

d4m5=d4mlf(@fes51,x,y,z,w,q);
disp(strcat('exaple Sec. 5.1 for m=',num2str(m)))
disp(d4m5)


%example of section 5.2
v=[1,2,3];
q=[1.1,2.2,3.3];
JFqv=d1mlf(@fes52,v,q); %Unlike the Fortran version, in Matlab 
                        %the dimension of f is not necessary.
                        %In the fortran version we decided to avoid one
                        %evaluation [dimf=size(f(q)), required by some 
                        %compilers], by passing the dimension of f 
                        %as argument.  

disp('exaple Sec. 5.2')
disp(JFqv)


%example of section 5.3
v=[0.1, 0.2];
q=[1.1+0.1.*1i, 2.2+0.1.*1i];
d3=d3mlf(@fes53,v,v,v,q);

disp('exaple Sec. 5.3')
disp(d3)

%computing the 4th derivative of the sine function (sin(q1*q2)) nested 
%500 times involve a considerable amout of operations (it takes 
%about 10 seconds)
%for heavy computations we recommend the Fortran version.
d4 = d4mlf(@fes53,v,v,v,v,q);
disp(d4)

q=[1.1,2.2];
q1p=[0.5, -2.7];
q2p=[-0.1, 0.7];

accel=d1mlf(@fes53,q2p,q) + d2mlf(@fes53,q1p,q1p,q);
disp(accel)


%example of section 5.4 (only using the dual number formulation)
%in many ocassions the partial derivatives can be avoided, using instead
%directional derivatives; which are computed much more efficiently
x0=[1.1, 2.2, 3.3, 4.4, 5.5];
EV=eye(5);
e1=EV(1,:); e2=EV(2,:); e3=EV(3,:); e4=EV(4,:); e5=EV(5,:);

disp('exaple Sec. 5.4')
dfx=d1mlf(@fes54,e1,x0);
disp(dfx)

dfzz=d2mlf(@fes54,e3,e3,x0);
disp(dfzz)

dfuwy=d3mlf(@fes54,e5,e4,e2,x0);
disp(dfuwy)

dfuzwx=d4mlf(@fes54,e5,e3,e4,e1,x0);
disp(dfuzwx)

%%%%%%%%%%%%%%%%%%%%%%%%%%
%functions
%example sec. 5.1
function retval = fes51(x)
    m = size(x,2);
	xcut1 = x(:,1:m-1);
	xcut2 = x(:,2:m); 
	num = -(xcut1.^2 + xcut2.^2 + 0.5*xcut1.*xcut2); 
	disc = xcut1.^2 + xcut2.^2 + 0.5.*xcut1.*xcut2;
	retval = -sum(exp(num./8).*cos(4*sqrt(disc)),2);  
end

%example sec. 5.2
function fr = fes52(q)
q1=q(1); q2=q(2); q3=q(3);
fr=[sin(q1.*q2.^2.*q3.^3), cos(q1.*q2.^2.*q3.^3)];
end

%example sec. 5.3
function fr = fes53(q)
q1=q(1); q2=q(2);

fr1=log(q1.*q2.^2); 
fr2=sqrt(q2)./q1;

fr3=sin(q1.*q2);
for k=1:500-1
    fr3=sin(fr3);
end
fr=cat(2,fr1,fr2,fr3);
end

%example sec. 5.4
function fr = fes54(q)
x=q(1); y=q(2); z=q(3); w=q(4); u=q(5);
fr=cos(x.*y./u).*z./w + 3.*sin(x.*u).*sin(y./u).*log(x.*y.*z./(u.*w));
end



