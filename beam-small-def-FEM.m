%{
1 Simply supported beam
2 E=1
3 I=1
4 L=1
5 equally divided elements n
6 Ramp loading with q=0 at x=0 ; q=qo at x=L;
7 The output will be displacement and slopes of all nodes with the exception of displacement of first and last nodes which are considered 0 in
boundary conditions 
%}


E= 1;
L=1;
r=0.1;
I=1;
n=5;    % Change the value of number of elements here
l=(L/n);
qo=1;

ke=[12,6*l,-12,6*l;6*l,4*l^2,-6*l,2*l^2;-12,-6*l,12,-6*l;6*l,2*l^2,-6*l,4*l^2;];
kg=zeros(2*n+2,2*n+2);

k=1;
while(k<=n)
p=kg(2*k-1,2*k-1);
q=ke(1,1);
kg(2*k-1,2*k-1)=(p+q);
p=kg(2*k-1,2*k);
q=ke(1,2);
kg(2*k-1,2*k)=(p+q);
p=kg(2*k-1,2*k+1);
q=ke(1,3);
kg(2*k-1,2*k+1)=(p+q);
p=kg(2*k-1,2*k+2);
q=ke(1,4);
kg(2*k-1,2*k+2)=(p+q);
p=kg(2*k,2*k-1);
q=ke(2,1);
kg(2*k,2*k-1)=(p+q);
p=kg(2*k,2*k);
q=ke(2,2);
kg(2*k,2*k)=(p+q);
p=kg(2*k,2*k+1);
q=ke(2,3);
kg(2*k,2*k+1)=(p+q);
p=kg(2*k,2*k+2);
q=ke(2,4);
kg(2*k,2*k+2)=(p+q);
p=kg(2*k+1,2*k-1);
q=ke(3,1);
kg(2*k+1,2*k-1)=(p+q);
p=kg(2*k+1,2*k);
q=ke(3,2);
kg(2*k+1,2*k)=(p+q);
p=kg(2*k+1,2*k+1);
q=ke(3,3);
kg(2*k+1,2*k+1)=(p+q);
p=kg(2*k+1,2*k+2);
q=ke(3,4);
kg(2*k+1,2*k+2)=(p+q);
p=kg(2*k+2,2*k-1);
q=ke(4,1);
kg(2*k+2,2*k-1)=(p+q);
p=kg(2*k+2,2*k);
q=ke(4,2);
kg(2*k+2,2*k)=(p+q);
p=kg(2*k+2,2*k+1);
q=ke(4,3);
kg(2*k+2,2*k+1)=(p+q);
p=kg(2*k+2,2*k+2);
q=ke(4,4);
kg(2*k+2,2*k+2)=(p+q);
k=(k+1);
end

Kg=E*I*(1/l^3)*kg;
Kg

fun1 = @(x,l) 1-3*x.^2/(l.^2)+2*x.^3/(l.^3);
q1= integral(@(x)fun1(x,l),0,l);

fun2 = @(x,l) x-2*x.^2/(l)+x.^3/(l.^2);
q2= integral(@(x)fun2(x,l),0,l);

fun3 = @(x,l) 3*x.^2/(l.^2)-2*x.^3/l.^3;
q3= integral(@(x)fun3(x,l),0,l);

fun4 = @(x,l) x.^3/(l^2)-1*x.^2/(l);
q4= integral(@(x)fun4(x,l),0,l);

fun5 = @(x,l) x-3*x.^3/(l.^2)+2*x.^4/(l.^3);
q5= integral(@(x)fun5(x,l),0,l);

fun6 = @(x,l) x.^2-2*x.^3/(l)+x.^4/(l.^2);
q6= integral(@(x)fun6(x,l),0,l);

fun7 = @(x,l) 3*x.^3/(l.^2)-2*x.^4/l.^3;
q7= integral(@(x)fun7(x,l),0,l);

fun8 = @(x,l) x.^4/(l^2)-1*x.^3/(l);
q8= integral(@(x)fun8(x,l),0,l);

Q1=[qo*q1; qo*q2; qo*q3; qo*q4;];
Q2=[qo*q5/L; qo*q6/L; qo*q7/L; qo*q8/L;];

Fext=zeros(2*n+2,1);

k=1;

while(k<=n)
    
s=(k-1)/n;
Q1a=s*Q1;
Q=Q1a+Q2;

p=Fext(2*k-1,1);
q=Q(1,1);  
Fext(2*k-1,1)=(p+q);


p=Fext(2*k,1);
q=Q(2,1);  
Fext(2*k,1)=(p+q);

p=Fext(2*k+1,1);
q=Q(3,1);  
Fext(2*k+1,1)=(p+q);

p=Fext(2*k+2,1);
q=Q(4,1);  
Fext(2*k+2,1)=(p+q);

k=(k+1);

end

Fext

KgBC=zeros(2*n,2*n);
i=2;
while(i<=2*n)
    j=2;
    while(j<=2*n)
        KgBC(i-1,j-1)=Kg(i,j);
        j=(j+1);
    end
    i=(i+1);
end   

i=2;
while(i<=2*n)
   KgBC(i-1,2*n)=Kg(i,2*n+2); 
    i=(i+1);
end

j=2;
while(j<=2*n+1)
    KgBC(2*n,j-1)=Kg(2*n+2,j);
    j=(j+1);
end

KgBC(2*n,2*n)=Kg(2*n+2,2*n+2);

KgBC


KgBCinv=inv(KgBC);
KgBCinv

FextBC=zeros(2*n,1);
i=2;
while(i<=2*n)
    FextBC(i-1,1)=Fext(i,1);
    i=(i+1);
end

FextBC(2*n,1)=Fext(2*n+2,1);

FextBC

u=KgBCinv*FextBC;
u
