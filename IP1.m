function [x,val,status]=IP1(f,A,b,Aeq,beq,lb,ub,M,e)
options = optimset('display','off');
bound=inf; % the initial bound is set to +ve infinity
[x0,val0]=linprog(f,A,b,Aeq,beq,lb,ub,[],options); 
[x,val,status,b]=rec(f,A,b,Aeq,beq,lb,ub,x0,val0,M,e,bound); % a recursive function that processes the BB tree 

function [xx,val,status,bb]=rec(f,A,b,Aeq,beq,lb,ub,x,v,M,e,bound) 
options = optimset('display','off');
% x is an initial solution and v is the corressponding objective function value

% solve the corresponding LP model with the integarily constraints removed
[x0,val0,status0]=linprog(f,A,b,Aeq,beq,lb,ub,[],options); 

% if the solution is not feasible or the value of the objective function is
% higher than the current bound return with the input intial solution
if status0<=0 || val0 > bound  
    xx=x; val=v; status=status0; bb=bound;
    return;
end

% if the integer-constraint variables turned to be integers within the
% input tolerance return
ind=find( abs(x0(M)-round(x0(M)))>e ); 
if isempty(ind)
    status=1;        
    if val0 < bound    % this solution is better than the current solution hence replace
        x0(M)=round(x0(M));
        xx=x0;        
        val=val0;
        bb=val0;
    else
        xx=x;  % return the input solution
        val=v;
        bb=bound;
    end
    return
end

% if we come here this means that the solution of the LP relaxation is
% feasible and gives a less value than the current bound but some of the
% integer-constraint variables are not integers. 
% Therefore we pick the first one that is not integer and form two LP problems
% and solve them recursively by calling the same function (branching)

% first LP problem with the added constraint that Xi <= floor(Xi) , i=ind(1)
br_var=M(ind(1));
br_value=x(br_var);
if isempty(A)
    [r c]=size(Aeq);
else
    [r c]=size(A);
end
A1=[A ; zeros(1,c)];
A1(end,br_var)=1;
b1=[b;floor(br_value)];

% second LP problem with the added constraint that Xi >= ceil(Xi) , i=ind(1)
A2=[A ;zeros(1,c)];
A2(end,br_var)=-1;
b2=[b; -ceil(br_value)];


% solve the first LP problem
[x1,val1,status1,bound1]=rec(f,A1,b1,Aeq,beq,lb,ub,x0,val0,M,e,bound);
status=status1;
if status1 >0 && bound1<bound % if the solution was successfull and gives a better bound
   xx=x1;
   val=val1;
   bound=bound1;
   bb=bound1;
else
    xx=x0;
    val=val0;
    bb=bound;
end
    
% solve the second LP problem
[x2,val2,status2,bound2]=rec(f,A2,b2,Aeq,beq,lb,ub,x0,val0,M,e,bound);

if status2 >0 && bound2<bound % if the solution was successfull and gives a better bound
    status=status2;
    xx=x2;
    val=val2;
    bb=bound2;
end
