clear all
close all
clc
for r = 5
%%Reading Files
name1 = ['nodeandcord', num2str(r),'.xlsx'];
name2 = ['elementnoderelation', num2str(r),'.xlsx'];
name4 = ['disp', num2str(r),'.xlsx'];
name5 = ['force', num2str(r),'.xlsx'];
ANTYP = input('Enter Analysis Type PLSTRESS: 21 ; PLSTRAIN: 22 - ');
TCK = input('Enter Thickness - ');
NODECOORD = xlsread(name1);
NNODE = length(NODECOORD); ANLYTYPE = NODECOORD(1:3,5);

ELMATNDRL = xlsread(name2);
NCA = ELMATNDRL(:,1:4);
[s,t] = size(NCA);
NELEM = s;
MTMAT = xlsread('mtrlmat1.xlsx');
[NMAT, hh] = size(MTMAT);

%% --------------------------------------------------------------------------------------------------- %%
                                        % START OF GEOMETRY HERE %
%-------------------------------------------------------------------------------------------------------%

% Coordinates
COORD(1:NNODE,1:3) = NODECOORD(1:NNODE,1:3); % Node no. 2nd Row X Cord, 3rd Y Cord

% Node Connection Array
for i = 1:NELEM
    XYCORDN1(i,1:2) = COORD(NCA(i,2),2:3);  % Here we take node no from NCA matrix ...            
    XYCORDN2(i,1:2) = COORD(NCA(i,3),2:3);  % and take x and y coord of respective ...            
    XYCORDN3(i,1:2) = COORD(NCA(i,4),2:3);  % node no for the given element no
    % row 2 is for y coord of XYCOORDELEMNj    % eg element no 1 has n1=1 n2=2 n3=4
    % row 1 is for x coord                     % so x coord of that is 0,0,2 shown in first ...
    % where j=1,2,3                            % of xcordelemn1,n2,and n3
end

% Making Mesh Diagram
for i = 1:NELEM
    XCORDN = [XYCORDN1(i,1); XYCORDN2(i,1); XYCORDN3(i,1); XYCORDN1(i,1)];
    YCORDN = [XYCORDN1(i,2); XYCORDN2(i,2); XYCORDN3(i,2); XYCORDN1(i,2)];
    line(XCORDN,YCORDN);
    hold on
end

%% --------------------------------------------------------------------------------------------------- %%
                                       % START OF ASSEMBLY HERE %
%-------------------------------------------------------------------------------------------------------%

GSTIFF = zeros(2*NNODE,2*NNODE);
for i = 1:NELEM
    AREAMAT(:,:,i) = [1, XYCORDN1(i,1), XYCORDN1(i,2); 1, XYCORDN2(i,1), XYCORDN2(i,2); ... 
              1, XYCORDN3(i,1), XYCORDN3(i,2)];     % Making 2DELTA Matrix
    TWODELTA(i) = det(AREAMAT(:,:,i));
       
    % BETAs
    BT1 = XYCORDN2(i,2) - XYCORDN3(i,2); BT2 = XYCORDN3(i,2) - XYCORDN1(i,2);
    BT3 = XYCORDN1(i,2) - XYCORDN2(i,2);
    % GAMMAs
    GM1 = XYCORDN3(i,1) - XYCORDN2(i,1); GM2 = XYCORDN1(i,1) - XYCORDN3(i,1);
    GM3 = XYCORDN2(i,1) - XYCORDN1(i,1);
    % B Matrix   
    B(:,:,i) = [BT1 0 BT2 0 BT3 0; 0 GM1 0 GM2 0 GM3; GM1 BT1 GM2 BT2 GM3 BT3]./(TWODELTA(i));
    % D Matrix
    NU = MTMAT(ELMATNDRL(i,5),3);       % Getting NU and E fro resprective material
    E = MTMAT(ELMATNDRL(i,5),2);
    if (ANTYP == 21)                    % Plain Stress D
    D(:,:,i) = (E/(1-(NU*NU))).*[1 NU 0; NU 1 0; 0 0 ((1-NU)/2)];
    elseif (ANTYP == 22)                % plain Strain D
    D(:,:,i) = (E/((1+NU)*(1-2*NU))).*[1-NU NU 0; NU 1-NU 0; 0 0 ((1-2*NU)/2)];
    end
    
    % K matrix elemental
    K(:,:,i) = (B(:,:,i)'*D(:,:,i)*B(:,:,i)).*(TCK).*(TWODELTA(i)/2);
    ESTIFF = K(:,:,i);
    
    % GSTIFF Filling
    RN = zeros(6,1);
    CN = zeros(6,1);
    RN(2*(1:3)-1) = 2*(NCA(i,2:4))-1;
    RN(2*(1:3)) = 2*(NCA(i,2:4));
    CN(2*(1:3)-1) = 2*(NCA(i,2:4))-1;
    CN(2*(1:3)) = 2*(NCA(i,2:4));
    
    for m = 1:6
        for n = 1:6
    GSTIFF(RN(m),CN(n)) = GSTIFF(RN(m),CN(n)) + ESTIFF(m,n);
        end
    end
end

%% --------------------------------------------------------------------------------------------------- %%
                                   % START OF BOUNDARY CONDITION HERE %
%-------------------------------------------------------------------------------------------------------%

DISP = xlsread(name4);
[mm,nn] = size(DISP); DSPELM = mm;
FORCE = xlsread(name5); [m,n] = size(FORCE);
FRELM = m;
FRC = zeros(2*NNODE,1);

for i = 1:FRELM
    if (FORCE(i,2) == 101)                  % 101 for Y direction only
        FRC(2*FORCE(i,1),1) = FORCE(i,4);
    elseif (FORCE(i,2) == 110)              % 110 for X direction only
        FRC(2*FORCE(i,1)-1,1) = FORCE(i,3);
    elseif (FORCE(i,2) == 111)              % 111 for both X and Y direction
        FRC(2*FORCE(i,1),1) = FORCE(i,4);
        FRC(2*FORCE(i,1)-1,1) = FORCE(i,3);
    end
end

for i = 1:DSPELM
    if (DISP(i,2) == 101)
        GSTIFF(2*DISP(i,1),2*DISP(i,1)) = 10^32;
        FRC(2*DISP(i,1),1) = DISP(i,4)*10^32;
    elseif (DISP(i,2) == 110)
        GSTIFF(2*DISP(i,1)-1,2*DISP(i,1)-1) = 10^32;
        FRC(2*DISP(i,1)-1,1) = DISP(i,3)*10^32;
    elseif (DISP(i,2) == 111)
        GSTIFF(2*DISP(i,1),2*DISP(i,1)) = 10^32;
        GSTIFF(2*DISP(i,1)-1,2*DISP(i,1)-1) = 10^32;
        FRC(2*DISP(i,1),1) = DISP(i,4)*10^32;
        FRC(2*DISP(i,1)-1,1) = DISP(i,3)*10^32;
    end
end
%% --------------------------------------------------------------------------------------------------- %%
                                % START OF POST PROCESSING CONDITION HERE %
%-------------------------------------------------------------------------------------------------------%

DSP(:,r) = linsolve(GSTIFF,FRC);

for i = 1:NELEM
    for j = 1:3
DSPELEM(2*j-1,i) = DSP(2*NCA(i,j+1)-1,r);
DSPELEM(2*j,i) = DSP(2*NCA(i,j+1),r);
    end
end

namevid = ['animdef', num2str(r)];
%V = VideoWriter(namevid,'Uncompressed AVI'); % initializing video making
%V.FrameRate = 10;                                % giving frame rate default is 30
%open(V)
ALPHA = linspace(0,1,100);

for j=1:100
    CORDNEW(1:NNODE,1) = 1:NNODE;
    CORDNEW(1:NNODE,2) = COORD(1:NNODE,2) + ALPHA(j).*DSP(2*(1:NNODE)-1,r);
    CORDNEW(1:NNODE,3) = COORD(1:NNODE,3) + ALPHA(j).*DSP(2*(1:NNODE),r);
    clf
    
    for i = 1:NELEM
        XYCORDN1NEW(i,1:2) = CORDNEW(NCA(i,2),2:3);           
        XYCORDN2NEW(i,1:2) = CORDNEW(NCA(i,3),2:3);            
        XYCORDN3NEW(i,1:2) = CORDNEW(NCA(i,4),2:3);   

    % Making Updated Mesh Diagram
        XCORDNEW = [XYCORDN1NEW(i,1); XYCORDN2NEW(i,1); XYCORDN3NEW(i,1); XYCORDN1NEW(i,1)];
        YCORDNEW = [XYCORDN1NEW(i,2); XYCORDN2NEW(i,2); XYCORDN3NEW(i,2); XYCORDN1NEW(i,2)];
        
     % Stress and Strain
        STRESS1(:,:,i) = D(:,:,i)*B(:,:,i)*DSPELEM(:,i);
        STRAIN1(:,:,i) = B(:,:,i)*DSPELEM(:,i);
        
        NU = MTMAT(ELMATNDRL(i,5),3);       % Getting NU and E fro resprective material
        E = MTMAT(ELMATNDRL(i,5),2);
        
        if (ANTYP == 21)                    % Plain Stress Epslion 33
            STRESS(:,:,i) = [STRESS1(1,1,i) STRESS1(3,1,i); STRESS1(3,1,i) STRESS1(2,1,i)];
            STRAIN(:,:,i) = [STRAIN1(1,1,i) STRAIN1(3,1,i) 0; STRAIN1(3,1,i) STRAIN1(2,1,i) 0; 0 0 0];
            STRAIN(3,3,i) = (-NU)/(1-NU)*(STRAIN(1,1,i)+STRAIN(2,2,i));
            MAXSTS(i,1) = max(STRESS(:,1,i));
            MAXSTN(i,1) = max(STRAIN(:,1,i));
            
            %Von mises stress
            PRNSTS(:,i) = eig(STRESS(:,:,i));
            VONSTS(i,1) = sqrt((PRNSTS(1,i))*(PRNSTS(1,i))+ ... 
            (PRNSTS(2,i))*(PRNSTS(2,i))-(PRNSTS(1,i))*(PRNSTS(2,i)));
        
        elseif (ANTYP == 22)                % plain Strain Sigma 33
            STRAIN(:,:,i) = [STRAIN1(1,1,i) STRAIN1(3,1,i); STRAIN1(3,1,i) STRAIN1(2,1,i)];
            STRESS(:,:,i) = [STRESS1(1,1,i) STRESS1(3,1,i) 0; STRESS1(3,1,i) STRESS1(2,1,i) 0; 0 0 0];
            STRESS(3,3,i) = (NU)*(STRESS(1,1,i)+STRESS(2,2,i));
            MAXSTS(i,1) = max(STRESS(:,1,i));
            MAXSTN(i,1) = max(STRAIN(:,1,i));
            
            PRNSTS(:,i) = eig(STRESS(:,:,i));
            VONSTS(i,1) = sqrt(0.5*((PRNSTS(1,i)-PRNSTS(2,i))^(2)+(PRNSTS(1,i)-PRNSTS(3,i))^(2)+ ...
                (PRNSTS(3,i)-PRNSTS(2,i))^(2)));
        end
        
        %Plots
        % Color Coding for different material
        if (ELMATNDRL(i,5) == 1)
            CLRDF = [1 0 0]; %Red
        elseif (ELMATNDRL(i,5) == 2)
            CLRDF = [0 0 1]; %Blue
        end
        
        %line(XCORDNEW,YCORDNEW,'Color','red','LineStyle','--');
        fill(XCORDNEW,YCORDNEW,VONSTS(i));
        pbaspect([1 1 1]);
        colorbar
        hold on
    end
    
%drawnow
pause(0.01)
%pic(j)=getframe(gcf);           % making plots into frame
end
%writeVideo(V,pic(1:100))        % writing frame in video
%close(gcf)                      %for making video of all pic gcf means current pic
%close(V)
end