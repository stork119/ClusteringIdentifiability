    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                                                           %
    %         Function includes system of ODEs describing       %
    %         NF-kB regulatory pathway.                         %
    %         Present molecules are coded as follows:           %
    %                                                           %  
    %         y(1)   IKKKa active                               %
    %         y(2)   IKKn   neutral                             %    
    %         y(3)   IKKa   active                              %
    %         y(4)   IKKi   inactive                            %
    %         y(5)   phospho-IkBa cytoplasmic                   %
    %         y(6)   (phospho-IkBa|NFkB) cytoplasmic            %  
    %         y(7)   NFkB  cytoplasmic                          %
    %         y(8)   NFkBn  nuclear                             %
    %         y(9)   A20                                        %
    %         y(10)  A20t                                       %
    %         y(11)  IkBa                                       %
    %         y(12)  IkBan                                      %
    %         y(13)  IkBat                                      %      
    %         y(14)  (IkBa|NFkB) cytoplasmic                    % 
    %         y(15)  (IkBan|NFkBn) nuclear                      %  
    %         y(16)  Active receptors                           %
    %         y(17)  A20 gene state                             %
    %         y(18)  IkBa gene state                            %
    %         y(19)  extracellular TNF                          % 
    %                                                           %
    %         y(20) Reporter gene state                         %
    %         y(21) Reporter transcript                         % 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 function ydot=ModelDclean(y, par, stim, t)
 % Ga,G,M,AN,ANa,ANR were taken from the input to the 'par' array
 %  y.length = 19, par.length = 56, t unused (but must be)

 par_cells = num2cell(par);
 [Ga,  G,    M,   AN,   ANa, ...
  ANR, NF0,  NF1, NF2,  M0,  ...
  M1,  M2,   k4,  ka20, AB,  ...
  kv,  q1,   q2,  c1,   c3,  ...
  c4,  c5,   k1,  k2,   k3,  ...
  a1,  a2,   a3,  c1a,  c5a, ...
  c6a, i1,   i1a, e1a,  e2a, ...
  dt,  tp,   KN,  KNN,  ka,  ...
  ki,  kb,   kf,  Tdeg, q1r, ...
  q2r, q2rr, c1r, c1rr, c3r, ...
  q1a, q2a,  c3a, c4a,  k5,  ...
  tp1 ] = par_cells{:};

 %###############################################################
 
 ydot(1)=ka*y(16)*(KN-y(1))* ka20/(ka20+y(9))-ki*y(1);                            %active IKKK kinase 
 ydot(2)=-y(1)^2*k1*y(2)+k5*(KNN-y(2)-y(3)-y(4));                                 %neutral IKK  (IKKn)
 ydot(3)=y(1)^2*k1*y(2)-k3*y(3)*(k2+y(9))/k2;                                     %free active IKK (IKKa)                                                                                     
 ydot(4)=k3*y(3)*(k2+y(9))/k2-k4*y(4);                                            %inactive IKK  (IKKi) 
 ydot(5)=a2*y(3)*y(11)-tp*y(5);                                                   %Phospo-IkBa cytoplasmic 
 ydot(6)=a3*y(3)*y(14)-tp1*y(6);                                                  %cytoplasmic (phospho-IkBa|NF-kB) 
 ydot(7)=c6a*y(14)-a1*y(7)*y(11)+tp1*y(6)-i1*y(7);                                %free cytoplasmic NFkB
 ydot(8)=i1*y(7)-a1*kv*y(12)*y(8);                                                %free nuclear NFkB
 ydot(9)=c4*y(10)-c5*y(9);                                                        %cytoplasmic A20
 ydot(10)=c1*y(17)-c3*y(10);                                                      %A20 transcript
 ydot(11)=-a2*y(3)*y(11)-a1*y(11)*y(7)+c4a*y(13)-c5a*y(11)-i1a*y(11)+e1a*y(12);   %free cytoplasmic IkBa
 ydot(12)=-a1*kv*y(12)*y(8)+i1a*y(11)-e1a*y(12);                                  %free nuclear IkBan
 ydot(13)=c1a*y(18)-c3a*y(13);                                                    %IkBa transcript
 ydot(14)=a1*y(11)*y(7)-c6a*y(14)-a3*y(3)*y(14)+e2a*y(15);                        %cytoplasmic (IkBa|NFkB) complex
 ydot(15)=a1*kv*y(12)*y(8)-e2a*y(15);                                             %nuclear (IkBa|NFkB) complex
 ydot(16)=kb*y(19)*(M-y(16))-kf*y(16);                                            %Active receptors
 ydot(17)=q1*y(8)*(AN-y(17))-q2*y(12)*y(17);                                      %A20 gene state
 ydot(18)=q1a*y(8)*(ANa-y(18))-q2a*y(12)*y(18);                                   %IkBa gene state
 ydot(19)=-Tdeg*y(19);                                                            %extracellular TNF
 %ydot(20)=q1r*y(8)*(ANR-y(20))-(q2rr+q2r*y(12))*y(20);                            % Reporter gene state 
 %ydot(21)=c1r*y(20)-c3r*y(21);                                                    % Reporter transcript  
   