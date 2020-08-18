close all
clc
format short g
roundn = @(t,n) round(t*10^n)./10^n;

n=24;            % Order of the graph

B= [[0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0], [1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [1, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0], [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0]];



A= reshape(B,[n,n]);
if A==A'
      disp('Matrix is Symmetric')
else
      disp('Not Symmetric')
end

q=eig(A);
x = roundn(q,4);
pi=sum(x>0);
ni=sum(x<0);
nullity=sum(x==0);
signature=abs(pi-ni);
rho=max(x);
E=sum(abs(eig(A)));
EE=sum(exp(eig(A)));

% The Laplacian matrix
l=size(A,1);
k=transpose(sum(A));
m=(sum(k))/2;
D = diag(k);
L=D-A;
q1=eig(L);
x1 = roundn(q1,4);
rho1=max(x1);
E1=sum(abs(eig(L)-((2*m)/n)));
EE1=sum(exp(eig(L)));

% The signless Laplacian matrix
Q=D+A;
q2=eig(Q);
x2 = roundn(q2,4);
rho2=max(x2);
E2=sum(abs(eig(Q)-((2*m)/n)));
EE2=sum(exp(eig(Q)));


% The extended adjacency matrix
A1=[];
for i=1:l
    for j=1:l
        if A(i,j)==1
            A(i,j)=(1/2)*((k(i)/k(j))+(k(j)/k(i)));
            A1=A;
        end
    end
end
q3=eig(A1);
x3 = roundn(q3,4);
rho3=max(x3);
E3=sum(abs(eig(A1)));


% The Randic matrix
R=[];
A= reshape(B,[n,n]);
for i=1:l
    for j=1:l
        if A(i,j)==1
            A(i,j)=1./sqrt(k(i)*k(j));
            R=A;
        end
    end
end
q4=eig(R);
x4 = roundn(q4,4);
rho4=max(x4);
E4=sum(abs(eig(R)));

% The general Randic matrix
% R1=[];
% syms a
% A= reshape(B,[n,n]);
% for i=1:l
%     for j=1:l
%         if A(i,j)==1
%             A(i,j)=(k(i)*k(j))^a;
%             R1=A;
%         end
%     end
% end
% Ra
% q5=eig(Ra);
% x5 = roundn(q5,4);
% rho5=max(x5);
% E5=sum(abs(eig(Ra)));

% The sum-connectivity matrix
S=[];
A= reshape(B,[n,n]);
for i=1:l
    for j=1:l
        if A(i,j)==1
            A(i,j)=1./sqrt(k(i)+k(j));
            S=A;
        end
    end
end
q5=eig(S);
x5 = roundn(q5,4);
rho5=max(x5);
E5=sum(abs(eig(S)));

% The ABC matrix
ABC=[];
A= reshape(B,[n,n]);
for i=1:l
    for j=1:l
        if A(i,j)==1
            A(i,j)=sqrt((k(i)+k(j)-2)/(k(i)*k(j)));
            ABC=A;
        end
    end
end
q6=eig(ABC);
x6 = roundn(q6,4);
rho6=max(x6);
E6=sum(abs(eig(ABC)));

% The GA matrix
GA=[];
A= reshape(B,[n,n]);
for i=1:l
    for j=1:l
        if A(i,j)==1
            A(i,j)=(2*sqrt(k(i)*k(j)))/((k(i)+k(j)));
            GA=A;
        end
    end
end
q7=eig(GA);
x7 = roundn(q7,4);
rho7=max(x7);
E7=sum(abs(eig(GA)));

% The AG matrix
AG=[];
A= reshape(B,[n,n]);
for i=1:l
    for j=1:l
        if A(i,j)==1
            A(i,j)=((k(i)+k(j)))/(2*sqrt(k(i)*k(j)));
            AG=A;
        end
    end
end
q8=eig(AG);
x8 = roundn(q8,4);
rho8=max(x8);
E8=sum(abs(eig(AG)));

% The first Zagreb matrix
Z1=[];
A= reshape(B,[n,n]);
for i=1:l
    for j=1:l
        if A(i,j)==1
            A(i,j)=(k(i)+k(j));
            Z1=A;
        end
    end
end
q9=eig(Z1);
x9 = roundn(q9,4);
rho9=max(x9);
E9=sum(abs(eig(Z1)));
EE9=sum(exp(eig(Z1)));

% The second Zagreb matrix
Z2=[];
A= reshape(B,[n,n]);
for i=1:l
    for j=1:l
        if A(i,j)==1
            A(i,j)=(k(i)*k(j));
            Z2=A;
        end
    end
end
q10=eig(Z2);
x10 = roundn(q10,4);
rho10=max(x10);
E10=sum(abs(eig(Z2)));
EE10=sum(exp(eig(Z2)));

% The harmonic matrix
H=[];
A= reshape(B,[n,n]);
for i=1:l
    for j=1:l
        if A(i,j)==1
            A(i,j)=2/(k(i)+k(j));
            H=A;
        end
    end
end
q11=eig(H);
x11 = roundn(q11,4);
rho11=max(x11);
E11=sum(abs(eig(H)));
EE11=sum(exp(eig(H)));



disp('Spectrum-related descriptors:')
fprintf('The A-SR is %4.4f\n',rho');
fprintf('The A-E is %4.4f\n',E');
fprintf('The A-EE is %4.4f\n',EE');
fprintf('The A-I+ is %d\n',pi');
fprintf('The A-I- is %d\n',ni');
fprintf('The A-I0 is %d\n',nullity');
% fprintf('The A-signature is %d\n',signature');
fprintf('The L-SR is %4.4f\n',rho1');
fprintf('The L-E is %4.4f\n',E1');
fprintf('The L-EE is %4.4f\n',EE1');
fprintf('The Q-SR is %4.4f\n',rho2');
fprintf('The Q-E is %4.4f\n',E2');
fprintf('The Q-EE is %4.4f\n',EE2');
fprintf('The A*-SR is %4.4f\n',rho3');
fprintf('The A*-E is %4.4f\n',E3');
fprintf('The R-SR is %4.4f\n',rho4');
fprintf('The R-E is %4.4f\n',E4');
fprintf('The S-SR is %4.4f\n',rho5');
fprintf('The S-E is %4.4f\n',E5');
fprintf('The ABC-SR is %4.4f\n',rho6');
fprintf('The ABC-E is %4.4f\n',E6');
fprintf('The GA-SR is %4.4f\n',rho7');
fprintf('The GA-E is %4.4f\n',E7');
fprintf('The AG-SR is %4.4f\n',rho8');
fprintf('The AG-E is %4.4f\n',E8');
fprintf('The Z1-SR is %4.4f\n',rho9');
fprintf('The Z1-E is %4.4f\n',E9');
fprintf('The Z1-EE is %4.4f\n',EE9');
fprintf('The Z2-SR is %4.4f\n',rho10');
fprintf('The Z2-E is %4.4f\n',E10');
fprintf('The Z2-EE is %4.4f\n',EE10');
fprintf('The H-SR is %4.4f\n',rho11');
fprintf('The H-E is %4.4f\n',E11');
fprintf('The H-EE is %4.4f\n',EE11');